# -*- coding:utf-8 -*-

import tables
import numpy as np


N_SIMU = 'n_simu'


def init_data_h5(filename):
    """Initialize a data HDF5 file"""
    with tables.openFile(filename, 'w') as f:
        f.createArray('/', N_SIMU, np.array(0))


def new_simu(filename, data):
    """Put the simulation data into the HDF5 file"""
    with tables.openFile(filename, 'a') as f:
        n_simu = f.getNode('/', N_SIMU).read()
        # parse data and put them in a new group
        simu_group = f.createGroup('/', 'simu' + str(n_simu))
        # TODO change value of n_simu
