#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from network_input import NetworkInput
from mitral_cells import MitralCells
from synapse import Synapse
from granule_cells import GranuleCells
from brian import *
import parameters as ps

import pickle
from multiprocessing import Pool

psmt = ps.Mitral()
psgr = ps.Granule()

def run((dt, trial)):
    print 'trial ', trial, 'dt =', dt
    
    # Clears Brian's cache and reinit simulation clock
    clear(all=True)
    defaultclock.reinit()

    # ## Network Input
    netin = NetworkInput()
    netin.set_eqs_model()
    netin.get_eqs_model()

    # ## Synapses
    synexc = Synapse(synapse_type='exc')
    synexc.set_eqs_model()
    synexc.get_eqs_model()

    syninhib = Synapse(synapse_type='inhib')
    syninhib.set_eqs_model()
    syninhib.get_eqs_model()

    # ## Mitral Cells
    N_mitral = 100
    mt = MitralCells()

    mt_supp_eqs = {'var': ['-I_input', '-I_syn'],
                   'eqs': [netin.get_eqs_model(), synexc.get_eqs_model()]}
    mt.add_eqs(supp_eqs=mt_supp_eqs)

    mt.make_pop(N_mitral)
    mt.pop.V = psmt.E_L

    # ## Granule Cells
    N_granule = 1
    gr = GranuleCells()

    gr_supp_eqs = {'var': ['-I_syn'],
                   'eqs': [syninhib.get_eqs_model()]}
    gr.add_eqs(supp_eqs=gr_supp_eqs)

    gr.make_pop(N_granule)
    gr.pop.V_D = psgr.E_L
    gr.pop.V_S = psgr.E_L

    # ## Connecting mitral and granule cells
    connections = ones((N_mitral,N_granule))
    @network_operation(when='start')
    def graded_synapse():
        #connections
        mt.pop.state('T')[:]=0.
        mt.pop.state('T')[mt.pop.get_refractory_indices()]=1.
        gr.pop.s_syn = dot(mt.pop.s, connections)
        mt.pop.s_syn = dot(gr.pop.s, transpose(connections))

    @network_operation(when='after_groups')
    def keep_reset():
        mt.pop.state('V')[mt.pop.get_refractory_indices()] = psmt.V_r

    # ## Simulation
    defaultclock.dt = dt
    t_simu = 2000*msecond

    monit_mitral = {}
    monit_granule = {}
    rec_mitral = [0, N_mitral/2, N_mitral-1]
    for pname in ('s', 's_syn', 'V', 'T', 'g_Ein'):
        monit_mitral[pname] = StateMonitor(mt.pop, pname, record=rec_mitral)
    for pname in ('V_D', 's_syn', 's', 'T'):
        monit_granule[pname] = StateMonitor(gr.pop, pname, record=True)
    monit_mitral['spikes'] = SpikeMonitor(mt.pop, record=True)

    nn = Network(mt.pop, gr.pop, [mgr for mgr in monit_granule.values()], [mmt for mmt in monit_mitral.values()],
                 graded_synapse, keep_reset)

    # nn.run(t_simu, report='text')
    nn.run(t_simu, report='text')

    # ## Results
    figure()
    plot(monit_mitral['s'].times/msecond, monit_mitral['s'][0], label="s mitral")
    plot(monit_mitral['s_syn'].times/msecond, monit_mitral['s_syn'][0], label="s_syn mitral")
    plot(monit_granule['s_syn'].times/msecond, monit_granule['s_syn'][0], label="s_syn granule")
    plot(monit_granule['s'].times/msecond, monit_granule['s'][0], label="s granule")
    legend()
    xlabel('time (ms)')
    ylabel('s mitral & s_syn granule & s granule')
    savefig('s_granule_dt'+str(dt)+'_trial'+str(trial)+'.pdf')

    raster_plot(monit_mitral['spikes'], newfigure=True)
    savefig('rasterplot_dt'+str(dt)+'_trial'+str(trial)+'.pdf')

    pickle.dump(monit_mitral['spikes'].spiketimes, open('spiketimes_dt'+str(dt)+'_trial'+str(trial)+'.txt', 'wb'))

if __name__ == '__main__':
    pool = Pool(processes=10)
    dts = [0.03*ms, 0.02*ms, 0.01*ms, 0.005*ms, 0.003*ms, 0.001*ms]
    inputs = []
    for dt in dts:
        for trial in xrange(3):
            inputs += [(dt, trial)]
    pool.map(run, inputs)
