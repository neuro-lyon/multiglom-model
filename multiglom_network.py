#!/usr/bin/env python2
# -*- coding:utf-8 -*-

"""
Multi-glomerular Network Simulation
===================================

This is a script that runs a simulation of a small OB model.
It consists of an input signal and some OB columns. A column is a set of one
glomerule, its connected mitral cells, and the granule cells.

Script Overview
---------------
1. Import the model and necessary stuff (like numpy)
2. Get the set of parameters (see inside the paramsets directory)
3. Initialize the different cell populations:
  - glomeruli
  - (synapses between granule and mitral)
  - mitral cells
  - granule cells
4. Connects the different cell populations:
  - glomeruli and mitral cells
  - mitral cells and granule cells
5. Set some monitors on the simulation
6. Run the simulation
7. Plots

"""

from brian import *
import numpy as np
from scipy.fftpack import fft, fftfreq

# One must import the 'parameters' module before any other modem.* import.
import model.parameters as ps

from model.glomerule import Glomerule
from model.mitral_cells import MitralCells
from model.synapse import Synapse
from model.granule_cells import GranuleCells


# Reset old stuff from Brian memory
clear(erase=True, all=True)
defaultclock.reinit()


"""
Parameters
----------
Get the parameter values from the `ps` module, which in turn gets the values
from the file specified in parameters.py.

Set some aliases for the different cell population sizes.
Also check that there is an even number of cells for each column.

Finally set some simulation parameters.

"""
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


"""
Population Initialization
-------------------------
1. glomeruli
*. synapses between granule and mitral cells
3. mitral cells
4. granule cells

"""
# Glomeruli
glom = Glomerule()
glom.add_eqs(oscillating=False)
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

# Granule Cells
gr = GranuleCells()
gr_supp_eqs = {'var': ['-I_syn'],
               'eqs': [syninhib.get_eqs_model()]}
gr.add_eqs(supp_eqs=gr_supp_eqs)

gr.make_pop(N_granule)
gr.pop.V_D = psgr.E_L
gr.pop.V_S = psgr.E_L


"""
Connecting Populations
----------------------
1. Glomeruli and mitral cells 
2. Mitral cells and granule cells

"""
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


"""
Simulation Monitoring
---------------------
Define a dict for each population. Then add the variable to monitor.

"""
# Simulation monitors
monit_glom = {}
monit_mt   = {}
monit_gr   = {}

recn = [0, N_mitral/2, N_mitral-1]
glom_pm = ('g')
mt_pm   = ('s', 's_syn', 'V')
gr_pm   = ('V_D', 's_syn', 's')
for pname in glom_pm:
    monit_glom[pname] = StateMonitor(glom.pop, pname, record=recn)
for pname in mt_pm:
    monit_mt[pname] = StateMonitor(mt.pop, pname, record=recn)
monit_mt['spikes'] = SpikeMonitor(mt.pop, record=True)
for pname in gr_pm:
    monit_gr[pname] = StateMonitor(gr.pop, pname, record=True)


"""
Running Simulation
------------------
Create Network object and put everything simulation related in it.
Then run this network.

"""
# Gathering simulation objects
netw = Network(glom.pop, mt.pop, gr.pop, mt_input, graded_synapse, keep_reset,
               [m for m in monit_glom.values()],
               [m for m in monit_mt.values()],
               [m for m in monit_gr.values()])

# Simulation run
netw.run(simu_length, report="text")


"""
Plotting
--------
Plot monitored variable and a scatter plot.

"""
# Raster plot
raster_plot(monit_mt['spikes'], newfigure=True)

# Membrane potentials
figure()
sub_v_mt = subplot(2, 1, 1)
for neur in recn:
    sub_v_mt.plot(monit_mt['V'].times/msecond,
                  monit_mt['V'][neur]/mvolt, label="mitral #" + str(neur))
sub_v_mt.legend()
sub_v_mt.set_xlabel('Time (ms)')
sub_v_mt.set_ylabel('Membrane potential of mitral : V (mvolt)')

sub_vd_gr = subplot(2, 1, 2, sharex=sub_v_mt)
for gran in xrange(N_granule):
    sub_vd_gr.plot(monit_gr['V_D'].times/msecond,
                   monit_gr['V_D'][gran]/mvolt, label="granule #" + str(gran))
sub_vd_gr.legend()
sub_vd_gr.set_xlabel('Time (ms)')
sub_vd_gr.set_ylabel('Membrane potential of granule : V (mvolt)')

# s and s_syn from granule and mitral cells
# also add an FFT on `s granule` to easily see the population frequency
figure()
sub_s = subplot(1, 2, 1)
sub_s.plot(monit_mt['s'].times/msecond,
         monit_mt['s'][0], label="s mitral")
sub_s.plot(monit_mt['s_syn'].times/msecond,
         monit_mt['s_syn'][0], label="s_syn mitral")
sub_s.plot(monit_gr['s_syn'].times/msecond,
         monit_gr['s_syn'][0], label="s_syn granule")
sub_s.plot(monit_gr['s'].times/msecond,
         monit_gr['s'][0], label="s granule")
sub_s.legend()
sub_s.set_xlabel('time (ms)')
sub_s.set_ylabel('s mitral & s_syn granule & s granule')

sub_syncrho = subplot(1, 2, 2)
FFT_MAX_FREQ = 200
NTIMES = len(monit_gr['s'].times)
FREQS = fftfreq(NTIMES, pscommon.simu_dt)
FFT_MAX_FREQ_INDEX = next(f for f in xrange(len(FREQS)) if FREQS[f] > FFT_MAX_FREQ)
sub_syncrho.plot(FREQS[:FFT_MAX_FREQ_INDEX],
     abs(fft(monit_gr['s'][0]-(monit_gr['s'][0]).mean())[:FFT_MAX_FREQ_INDEX]))
sub_syncrho.set_xlabel("granule 's' frequency (Hz)")
sub_syncrho.set_ylabel('Power')

show()
