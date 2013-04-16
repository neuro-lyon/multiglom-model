# -*- coding:utf-8 -*-

import tables
import os


def init_data_h5(filename):
    """Initialize a data HDF5 file"""
    # Raise en error if the file already exists
    try:
        os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0644)
    except OSError, e:
        raise e
    # Else, continue by creating the file
    else:
        with tables.openFile(filename, 'w') as f:
            setattr(f.root._v_attrs, 'n_simu', 0)


def new_simu(filename, data):
    """Put the simulation data into the HDF5 file"""
    with tables.openFile(filename, 'a') as f:
        n_simu = getattr(f.root._v_attrs, 'n_simu')
        # parse data and put them in a new group
        simu_group = f.createGroup('/', 'simu' + str(n_simu))
        # TODO change value of n_simu
