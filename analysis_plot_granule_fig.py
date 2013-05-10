# -*- coding:utf-8 -*-

"""
Select a specific simulation from a HDF5 file and plot its granule figure.

"""

import numpy as np
import matplotlib.pyplot as plt
import tables

import h5manager as h5m
from plotting import granule_pop_figure
from analysis import fftmax

DB_FILENAME = "data/db30x30_two_glom_beta_new_ps_interco_strength0_1_interco_rate0_1.h5"
DB = tables.openFile(DB_FILENAME)

# Get all simulation data
ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
         ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
         ('paramset', 'arrays', 'times'),
         ('paramset', '_v_attrs', 'Common', 'simu_dt'),
         ('results', 's_granule'),
         ('results', 's_syn_self'),
         ('results', '_v_attrs'))
ALL_SIMU_ATTRS = h5m.get_all_attrs(DB, ATTRS)

SET_INTERCO_RATE = list(set(s[0] for s in ALL_SIMU_ATTRS))
SET_INTERCO_RATE.sort()
SET_INTERCO_STRENGTH = list(set(s[1] for s in ALL_SIMU_ATTRS))
SET_INTERCO_STRENGTH.sort()


# Some useful functions
def deq(a, b, delta):
    """Delta equality"""
    return abs(a - b) < delta

class DummyPkg:
    def __init__(self, values, times):
        self.values = values
        self.times = times

# Get the specific simulations
def get_interco(simu, interco_rate, interco_strength):
    return deq(simu[0], interco_rate, 0.02) and deq(simu[1], interco_strength, 0.02)

INTERCO_RATE = SET_INTERCO_RATE[25]
INTERCO_STRENGTH = SET_INTERCO_STRENGTH[4]
GOOD_SIMUS = filter(lambda x : get_interco(x, INTERCO_RATE, INTERCO_STRENGTH), ALL_SIMU_ATTRS)

REDO_FFTMAX = []
for simu in GOOD_SIMUS:
    gr_s = simu[4].read()
    gr_s_syn_self = simu[5].read()
    times = simu[2].read()
    dt = float(simu[3])
    granule_pop_figure(gr_s, gr_s_syn_self, times, dt)

    dummy = DummyPkg(gr_s, times)
    REDO_FFTMAX.append(fftmax(dummy, 2, dt))

print REDO_FFTMAX
plt.show()

# Finally close the DB
DB.close()
