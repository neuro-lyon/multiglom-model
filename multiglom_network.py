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
APARSER.add_argument('--full-ps', action='store_true')
ARGS = APARSER.parse_args()

# Set the parameters from the specified file BEFORE any model.* import
from brian import *
import model_utils as mutils
import model
mutils.set_model_ps(ARGS.psfile)

import numpy as np
from scipy.fftpack import fft, fftfreq
import analysis
import plotting
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
mtgr_connections = mutils.intrapop_connections(N_MITRAL, N_GRANULE, N_SUBPOP, N_MITRAL_PER_SUBPOP)

# Inter subpopulation connectivities
INTER_CONN_RATE = PSCOMMON['inter_conn_rate']
INTER_CONN_STRENGTH = PSCOMMON['inter_conn_strength']
mtgr_connections = mutils.interpop_connections(mtgr_connections, N_MITRAL, N_SUBPOP,
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
monit_glom = mutils.monit(glom.pop, glom_ps, rec_neurons)
monit_mt   = mutils.monit(mt.pop, mt_ps, rec_neurons, spikes=True)
monit_gr   = mutils.monit(gr.pop, gr_ps)


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

if ARGS.full_ps:
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
    plotting.raster_plot(monit_mt['spikes'], N_SUBPOP)

    # Membrane potentials
    plotting.memb_plot_figure(monit_mt, monit_gr, rec_neurons, N_GRANULE)

    # Granule synapses
    plotting.granule_figure(monit_gr, PSCOMMON)

    show()
