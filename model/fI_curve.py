#!/usr/bin/env python2
# -*- coding:utf-8 -*-

"""
Script to plot the f-I curve of the mitral cells.
"""

from sys import argv    
from brian import *
from mitral_cells import MitralCells

start_I = float(argv[1])
stop_I = float(argv[2])
N_mitral = float(argv[3])
period = float(argv[4])*second

def make_plot_fI(spiketimes, currents, period):
    freqs = []
    print 'stimes', spiketimes
    for neuron in spiketimes:
        print 'n',neuron,'nspikes',len(spiketimes[neuron])
        freqs += [len(spiketimes[neuron])/period]
    plot(array(currents)/(amp*meter**-2), array(freqs)/second**-1)
    xlabel('current I (amp*meter**-2)')
    ylabel('frequency f (Hz)')
    show()

if __name__ == '__main__':
    # Initialize an empty mitral cell population
    mt = MitralCells()
    # Add electrode current
    supp = {'var': ['+I_electrode'], 'eqs' : ['I_electrode : amp*meter**-2']}
    mt.add_eqs(supp)
    # Make the population
    mt.make_pop(N_mitral)
    # Set the current for each mitral
    currents = linspace(start_I, stop_I, N_mitral)
    mt.pop.I_electrode = currents
    # Set the monitor to record the spikes
    monitor_spikes = SpikeMonitor(mt.pop, record=True)
    # Make the network and run the simulation
    netw = Network(mt.pop,  monitor_spikes)
    netw.run(period, report='text')
    # Make the fI curve and plot the result
    raster_plot(monitor_spikes)
    figure()
    make_plot_fI(monitor_spikes.spiketimes, currents, period)
