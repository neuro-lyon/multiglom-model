#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import *
import numpy as np

from model.glomerule import Glomerule
from model.network_input import NetworkInput
from model.mitral_cells import MitralCells
from model.synapse import Synapse
from model.granule_cells import GranuleCells
import model.parameters as ps

# Parameters
psgm     = ps.Glomerule()
psmt     = ps.Mitral()
psgr     = ps.Granule()
pscommon = ps.Common()

N_mitral    = pscommon.N_mitral
N_glomerule = N_granule = N_subpop = pscommon.N_subpop

# check to have an even number of mitral in each sub-population
assert N_mitral % N_subpop == 0, \
       "N_mitral is not a multiple of the number of sub-populations N_subpop."
N_mitral_per_subpop = N_mitral/N_subpop

defaultclock.dt = pscommon.simu_dt
simu_length     = pscommon.simu_length

# Glomeruli
glom = Glomerule()
glom.add_eqs()
glom.make_pop(N_glomerule*N_mitral_per_subpop)

# Synapses (granule -- mitral)
synexc = Synapse(synapse_type='exc') # excitatory synapse
synexc.set_eqs_model()
synexc.get_eqs_model()

syninhib = Synapse(synapse_type='inhib') # inhibitory synapse
syninhib.set_eqs_model()
syninhib.get_eqs_model()

# Mitral cells
mt = MitralCells()
mt_supp_eqs =  {'var': ['- I_syn', '- g_input*V'],
                'eqs': [synexc.get_eqs_model(),
                        Equations("g_input : siemens*meter**-2")]}
mt.add_eqs(supp_eqs=mt_supp_eqs)

mt.make_pop(N_mitral)
mt.pop.V = (psmt.E_L - psmt.V_r)*np.random.random_sample(np.shape(mt.pop.V)) \
           + psmt.V_r

# Mitral sub-populations
#mtsubpop = []
#for subpop in xrange(N_subpop):
    #mtsubpop += [mt.pop.subgroup(N_mitral/N_subpop)]

# Granule Cells
gr = GranuleCells()
gr_supp_eqs = {'var': ['-I_syn'],
               'eqs': [syninhib.get_eqs_model()]}
gr.add_eqs(supp_eqs=gr_supp_eqs)

gr.make_pop(N_granule)
gr.pop.V_D = psgr.E_L
gr.pop.V_S = psgr.E_L

# Connecting mitral cells to glomeruli
glmt_connections = diag(ones(N_mitral))

# Glomeruli--Mitral interactions
@network_operation(when='start')
def mt_input():
    mt.pop.g_input = dot(glom.pop.g, glmt_connections)

# Connecting sub-population of mitral cells to granule cells
mtgr_connections = np.zeros((N_mitral, N_granule))
for i_subpop in xrange(N_subpop):
    start = i_subpop*N_mitral_per_subpop
    stop  = start + N_mitral_per_subpop
    mtgr_connections[start:stop, i_subpop] = 1.

# Mitral--Granule interactions
@network_operation(when='start')
def graded_synapse():
    mt.pop.state('T')[:] = 0.
    mt.pop.state('T')[mt.pop.get_refractory_indices()] = 1.
    gr.pop.s_syn = dot(mt.pop.s, mtgr_connections)
    mt.pop.s_syn = dot(gr.pop.s, transpose(mtgr_connections))

@network_operation(when='after_groups')
def keep_reset():
    mt.pop.state('V')[mt.pop.get_refractory_indices()] = psmt.V_r


# Simulation monitors
monit_glom = {}
monit_mt   = {}
monit_gr   = {}

recn = [0, N_mitral/2, N_mitral-1]
glom_pm = ('g')
mt_pm   = ('s', 's_syn', 'V', 'T')
gr_pm   = ('V_D', 's_syn', 's', 'T')
for pname in glom_pm:
    monit_glom[pname] = StateMonitor(glom.pop, pname, record=recn)
for pname in mt_pm:
    monit_mt[pname] = StateMonitor(mt.pop, pname, record=recn)
monit_mt['spikes'] = SpikeMonitor(mt.pop, record=True)
for pname in gr_pm:
    monit_gr[pname] = StateMonitor(gr.pop, pname, record=True)

# Gathering simulation objects
netw = Network(glom.pop, mt.pop, gr.pop, mt_input, graded_synapse, keep_reset,
               [m for m in monit_glom.values()],
               [m for m in monit_mt.values()],
               [m for m in monit_gr.values()])

# Simulation run and plots
netw.run(simu_length, report="text")

figure()
for n in recn:
    plot(monit_glom['g'].times/msecond,
         monit_glom['g'][n]/mvolt, label="glom neur #"+str(n))
legend()
xlabel('time (ms)')
ylabel('g of glomeruli (siemens*amp*meter**-2)')

figure()
for n in recn:
    plot(monit_mt['V'].times/msecond,
         monit_mt['V'][n]/mvolt, label="mitral #"+str(n))
legend()
xlabel('time (ms)')
ylabel('membrane potential of mitral : V (mvolt)')

figure()
plot(monit_mt['s'].times/msecond,
     monit_mt['s'][0], label="s mitral")
plot(monit_mt['s_syn'].times/msecond,
     monit_mt['s_syn'][0], label="s_syn mitral")
plot(monit_gr['s_syn'].times/msecond,
     monit_gr['s_syn'][0], label="s_syn granule")
plot(monit_gr['s'].times/msecond,
     monit_gr['s'][0], label="s granule")
legend()
xlabel('time (ms)')
ylabel('s mitral & s_syn granule & s granule')

print 'Spikes of', N_mitral, ' mitral:', monit_mt['spikes'].nspikes
raster_plot(monit_mt['spikes'], newfigure=True)

show()
