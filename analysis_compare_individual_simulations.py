# -*- coding:utf-8 -*-

"""
Select specific simulations from a HDF5 file and plot the granule figures and
scatter plot to compare the simulations.

"""

import matplotlib.pyplot as plt
import tables
from numpy import allclose

import h5manager as h5m
from plotting import granule_pop_figure
from analysis import fftmax

DB_FILENAME = "data/db30x30_beta_homeostasis.h5"
DB = tables.openFile(DB_FILENAME)

# Get all simulation data
ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
         ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
         ('paramset', 'arrays', 'times'),
         ('paramset', '_v_attrs', 'Common', 'simu_dt'),
         ('results', 's_granule'),
         ('results', 's_syn_self'),
         ('results', '_v_attrs'),
         ('results', 'spikes_it'),
)
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

        gr_s = simu[4].read()
        gr_s_syn_self = simu[5].read()
        times = simu[2].read()
        dt = float(simu[3])
        granule_pop_figure(gr_s, gr_s_syn_self, times, dt)

        signal = SignalRepack(gr_s, times)
        REDO_FFTMAX.append(fftmax(signal, 2, dt))

        spikes_it = simu[7].read()
        plt.figure()
        plt.title("rate: "+str(rate)+", strength: "+str(strength))
        plt.plot(spikes_it[1], spikes_it[0], ' .', mew=0)

        print 'rate:', rate, 'strength:', strength, REDO_FFTMAX
plt.show()

# Finally close the DB
DB.close()
