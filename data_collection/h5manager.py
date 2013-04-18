# -*- coding:utf-8 -*-

import tables
import os


def init_data_h5(filename):
    """Initialize a data HDF5 file"""
    filename = filename_to_h5(filename)
    if not file_exists(filename):
        f = tables.openFile(filename, 'w')
        f.close()


def write_simu_data(filename, info, paramset, results):
    """Create a HDF5 file for new simulation"""
    with tables.openFile(filename, 'a') as f:
        # Put info into the HDF5 root
        for attr in info:
            setattr(f.root._v_attrs, attr, info[attr])
        # Put the parameter set into the file
        pset = f.createGroup('/', 'paramset', title="Parameter set")
        for attr in paramset:
            setattr(pset._v_attrs, attr, paramset[attr])
        # Put the data results into the file
        res = f.createGroup('/', 'results', title="Simulation results")
        for attr in results:
            f.createArray(res, attr, results[attr][0], title=results[attr][1])


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


def filename_to_h5(filename):
    """Add extension .h5 to filename"""
    return filename if filename.endswith('.h5') else filename + '.h5'
