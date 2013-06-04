# -*- coding:utf-8 -*-

"""
Select specific simulations from a HDF5 file and plot the granule figures and
scatter plot to compare the simulations.

"""

import matplotlib.pyplot as plt
import tables
from numpy import allclose

import h5manager as h5m
from plotting import granule_pop_figure, raster_plot
from analysis import fftmax
from arg_parsers import ANACOMP_PARSER

ARGS = ANACOMP_PARSER.parse_args()

DB_FILENAME = ARGS.data_file
DB = tables.openFile(DB_FILENAME)

PLOT_MEMB_POT = ARGS.plot_mp

# Get all simulation data
ATTRS = [('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
         ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
         ('paramset', 'arrays', 'times'),
         ('paramset', '_v_attrs', 'Common', 'simu_dt'),
         ('results', 's_granule'),
         ('results', 's_syn_self'),
         ('results', '_v_attrs'),
         ('results', 'spikes_it'),
         ('paramset', 'arrays', 'mtgr_connections'),
]
if PLOT_MEMB_POT:
    ATTRS.append(('results', 'mean_memb_pot'))

ALL_SIMU_ATTRS = h5m.get_all_attrs(DB, ATTRS)
ALL_SIMU_ATTRS.sort()

SET_INTERCO_RATE = list(set(s[0] for s in ALL_SIMU_ATTRS))
SET_INTERCO_RATE.sort()
SET_INTERCO_STRENGTH = list(set(s[1] for s in ALL_SIMU_ATTRS))
SET_INTERCO_STRENGTH.sort()


# Some useful functions
def deq(a, b, delta):
    """Delta equality"""
    return abs(a - b) < delta

class SignalRepack:
    def __init__(self, values, times):
        self.values = values
        self.times = times

# Get the specific simulations
def get_interco(simu, interco_rate, interco_strength):
    return deq(simu[0], interco_rate, 0.02) and deq(simu[1], interco_strength, 0.02)

SELECTED_RATES = [20, 25]
SELECTED_STRENGTH = [6]
REDO_FFTMAX = []
for rate in SELECTED_RATES:
    for strength in SELECTED_STRENGTH:
        interco_rate = SET_INTERCO_RATE[rate]
        interco_strength = SET_INTERCO_STRENGTH[strength]
        filtered_simu = filter(lambda x: allclose(x[0], interco_rate) \
                                         and allclose(x[1], interco_strength),
                               ALL_SIMU_ATTRS)
        if len(filtered_simu) > 1:
            print "Warning: more than one simulation satisfy your criteria,",
            print "selecting the first one."
        simu = filtered_simu[0]

        # Granule plot
        gr_s = simu[4].read()
        gr_s_syn_self = simu[5].read()
        times = simu[2].read()
        dt = float(simu[3])
        mtgr_connections = simu[8].read()
        granule_pop_figure(gr_s, gr_s_syn_self, times, dt)

        # Raster plot
        spikes_it = simu[7].read()
        raster_plot(spikes_it[0], spikes_it[1], mtgr_connections)

        # Membrane potential
        if PLOT_MEMB_POT:
            memb_potentials = simu[9].read()
            labels = ("Mean. non-interco.", "Mean. interco")
            plt.figure()
            for ind_mp, memb_pot in enumerate(memb_potentials):
                label = labels[ind_mp % 2]
                label += " glom #" + str(ind_mp/2)
                plt.plot(times, memb_pot, label=label)
                plt.legend()
            plt.xlabel("Time (s)")
            plt.ylabel("Membrane potential (V)")

        # FFT max peak
        signal = SignalRepack(gr_s, times)
        REDO_FFTMAX.append(fftmax(signal, 2, dt))

        print 'rate:', rate, 'strength:', strength, REDO_FFTMAX
plt.show()

# Finally close the DB
DB.close()
