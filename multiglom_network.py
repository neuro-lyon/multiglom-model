#!/usr/bin/env python2.7
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
7. Output simulation information and indexes.
8. Plots

"""

# Argument parsing
import argparse

APARSER = argparse.ArgumentParser(description="Run a multi-glomerular simulation.")
APARSER.add_argument('psfile')
APARSER.add_argument('--no-plot', action='store_true')
APARSER.add_argument('--no-indexes', action='store_true')
APARSER.add_argument('--no-full-ps', action='store_true')
ARGS = APARSER.parse_args()

# Set the parameters from the specified file BEFORE any model.* import
from brian import *
from model_utils import set_model_ps, intrapop_connections, interpop_connections, monit

import model

set_model_ps(ARGS.psfile)

import numpy as np
from scipy.fftpack import fft, fftfreq
import analysis
from utils import print_dict

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
PSGM     = model.PARAMETERS['Glomerule']
PSMT     = model.PARAMETERS['Mitral']
PSGR     = model.PARAMETERS['Granule']
PSCOMMON = model.PARAMETERS['Common']

N_MITRAL    = PSCOMMON['N_mitral']
N_GLOMERULI = N_GRANULE = N_SUBPOP = PSCOMMON['N_subpop']

# check to have an even number of mitral in each sub-population
assert N_MITRAL % N_SUBPOP == 0, \
       "N_mitral is not a multiple of the number of sub-populations N_subpop."
N_MITRAL_PER_SUBPOP = N_MITRAL/N_SUBPOP

defaultclock.dt = PSCOMMON['simu_dt']
SIMU_LENGTH     = PSCOMMON['simu_length']


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
glom.make_pop(N_GLOMERULI*N_MITRAL_PER_SUBPOP)

# Synapses (granule -- mitral)
synexc = Synapse(synapse_type='exc') # excitatory synapse
synexc.set_eqs_model()

syninhib = Synapse(synapse_type='inhib') # inhibitory synapse
syninhib.set_eqs_model()

# Mitral cells
mt = MitralCells()
mt_supp_eqs =  {'var': ['- I_syn', '- g_input*V'],
                'eqs': [synexc.get_eqs_model(),
                        Equations("g_input : siemens*meter**-2")]}
mt.add_eqs(supp_eqs=mt_supp_eqs)
mt.make_pop(N_MITRAL)
mt.pop.V = (PSMT['E_L'] - PSMT['V_r'])*np.random.random_sample(np.shape(mt.pop.V)) \
           + PSMT['V_r']

# Granule Cells
gr = GranuleCells()
gr_supp_eqs = {'var': ['-I_syn'],
               'eqs': [syninhib.get_eqs_model()]}
gr.add_eqs(supp_eqs=gr_supp_eqs)
gr.make_pop(N_GRANULE)
gr.pop.V_D = PSGR['E_L']
gr.pop.V_S = PSGR['E_L']


"""
Connecting Populations
----------------------
1. Glomeruli and mitral cells 
2. Mitral cells and granule cells

"""
# Connecting mitral cells to glomeruli
glmt_connections = diag(ones(N_MITRAL))

# Glomeruli--Mitral interactions
@network_operation(when='start')
def mt_input():
    mt.pop.g_input = dot(glom.pop.g, glmt_connections)

# Connecting sub-population of mitral cells to granule cells
mtgr_connections = intrapop_connections(N_MITRAL, N_GRANULE, N_SUBPOP, N_MITRAL_PER_SUBPOP)

# Inter subpopulation connectivities
INTER_CONN_RATE = PSCOMMON['inter_conn_rate']
INTER_CONN_STRENGTH = PSCOMMON['inter_conn_strength']
mtgr_connections = interpop_connections(mtgr_connections, N_MITRAL, N_SUBPOP,
                                        N_MITRAL_PER_SUBPOP, INTER_CONN_RATE, INTER_CONN_STRENGTH)
# Mitral--Granule interactions
@network_operation(when='start')
def graded_synapse():
    mt.pop.state('T')[:] = 0.
    mt.pop.state('T')[mt.pop.get_refractory_indices()] = 1.
    gr.pop.s_syn = dot(mt.pop.s, mtgr_connections)
    mt.pop.s_syn = dot(gr.pop.s, transpose(mtgr_connections))

@network_operation(when='after_groups')
def keep_reset():
    mt.pop.state('V')[mt.pop.get_refractory_indices()] = PSMT['V_r']


"""
Simulation Monitoring
---------------------
Monitor state variables for the different populations.

"""
rec_neurons = [0, N_MITRAL/2, N_MITRAL - 1]
glom_ps = ('g')
mt_ps   = ('s', 's_syn', 'V')
gr_ps   = ('V_D', 's_syn', 's')

# Simulation monitors
monit_glom = monit(glom.pop, glom_ps, rec_neurons)
monit_mt   = monit(mt.pop, mt_ps, rec_neurons, spikes=True)
monit_gr   = monit(gr.pop, gr_ps)


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
netw.run(SIMU_LENGTH, report="text")


"""
Information Output
------------------

"""
print '\nParameters: using', ARGS.psfile

print 'Populations:', N_SUBPOP, 'glomerular columns;',
print N_MITRAL, 'mitral cells;', N_GRANULE, 'granule cells.'

print 'Times:', SIMU_LENGTH, 'of simulation; dt =', defaultclock.dt, '.'

if not ARGS.no_full_ps:
    print 'Full set of parameters:'
    print_dict(model.PARAMETERS)

if not ARGS.no_indexes:
    sts_index = analysis.sts(monit_gr['s_syn'], monit_mt['spikes'])
    mps_index = analysis.mps(monit_mt['V'])
    print 'Indexes: STS =', sts_index, '; MPS =', mps_index, '.'


"""
Plotting
--------
Plot monitored variables and a scatter plot.

"""
if not ARGS.no_plot:
    # Raster plot
    raster_plot(monit_mt['spikes'], newfigure=True)

    # Membrane potentials
    figure()
    sub_v_mt = subplot(2, 1, 1)
    for neur in rec_neurons:
        sub_v_mt.plot(monit_mt['V'].times/msecond,
                      monit_mt['V'][neur]/mvolt)
    sub_v_mt.set_xlabel('Time (ms)')
    sub_v_mt.set_ylabel('Membrane potential of mitral : V (mvolt)')

    sub_vd_gr = subplot(2, 1, 2, sharex=sub_v_mt)
    for gran in xrange(N_GRANULE):
        sub_vd_gr.plot(monit_gr['V_D'].times/msecond,
                       monit_gr['V_D'][gran]/mvolt, label="granule #" + str(gran))
    sub_vd_gr.legend()
    sub_vd_gr.set_xlabel('Time (ms)')
    sub_vd_gr.set_ylabel('Membrane potential of granule : V (mvolt)')

    # s and s_syn from granule and mitral cells
    # also add an FFT on `s granule` to easily see the population frequency
    for gr in xrange(N_GRANULE):
        figure()
        sub_s = subplot(1, 2, 1)
        sub_s.plot(monit_gr['s_syn'].times/msecond,
                 monit_gr['s_syn'][gr], label="s_syn granule #"+str(gr))
        sub_s.plot(monit_gr['s'].times/msecond,
                 monit_gr['s'][gr], label="s granule #"+str(gr))
        sub_s.legend()
        sub_s.set_xlabel('time (ms)')
        sub_s.set_ylabel('s mitral & s_syn granule & s granule #'+str(gr))

        sub_syncrho = subplot(1, 2, 2)
        fft_max_freq = 200
        ntimes = len(monit_gr['s'].times)
        freqs = fftfreq(ntimes, PSCOMMON['simu_dt'])
        fft_max_freq_index = next(f for f in xrange(len(freqs)) if freqs[f] > fft_max_freq)

        fft_sig = abs(fft(monit_gr['s'][gr]-(monit_gr['s'][gr]).mean())[:fft_max_freq_index])
        ind_max_freq = argmax(fft_sig)
        print 'MAX Freq FFT :', freqs[ind_max_freq]

        sub_syncrho.plot(freqs[:fft_max_freq_index], fft_sig)
        sub_syncrho.set_xlabel("granule #"+str(gr)+" 's' frequency (Hz)")
        sub_syncrho.set_ylabel('Power')

    show()
