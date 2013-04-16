# -*- coding:utf-8 -*-

import tables
import os


def init_data_h5(filename):
    """Initialize a data HDF5 file"""
    if not file_exists(filename):
        with tables.openFile(filename, 'w') as f:
            setattr(f.root._v_attrs, 'n_simu', 0)


def file_exists(filename):
    """Check if a file `filename` exists."""
    file_exists = False
    try:
        fd = os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0644)
        os.close(fd)
        os.remove(filename)
    except OSError, e:
        file_exists = True
        raise e
    return file_exists


def new_simu(filename, data):
    """Put the simulation data into the HDF5 file"""
    with tables.openFile(filename, 'a') as f:
        n_simu = getattr(f.root._v_attrs, 'n_simu')
        # parse data and put them in a new group
        simu_group = f.createGroup('/', 'simu' + str(n_simu))
        # TODO change value of n_simu
