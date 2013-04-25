# -*- coding:utf-8 -*-

import tables
import os
from os import path
from utils import listdir_filter


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
        for attr in results['data']:
            f.createArray(res, attr, results['data'][attr][0],
                          title=results['data'][attr][1])
        for attr in results['indexes']:
            setattr(res._v_attrs, attr, results['indexes'][attr])


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


def collect_h5_to_db(dirpath, dbpath):
    """Collects all HDF5 in dirpath and put them into the big HDF5 DB"""
    h5files = listdir_filter(dirpath, lambda fname: fname[-3:] == '.h5')
    h5files.remove(dbpath)  # In case db is in the same directory
    old_simus_id = get_all_simus_id(dbpath)

    # Open the db to put the new simus in
    with tables.openFile(dbpath, 'a') as db:
        for h5file in h5files:
            # Browse all simulation files in dirpath
            with tables.openFile(h5file, 'r') as f:
                file_time = f.root._v_attrs['time']
                file_uuid = f.root._v_attrs['uuid']

                # Add the simulation if not in DB
                if (file_time, file_uuid) not in old_simus_id:
                    simu_id = file_time + '__' + file_uuid
                    simu_id = simu_id.replace(':', '_')
                    simu_id = simu_id.replace('.', '_')
                    simu_id = simu_id.replace('-', '_')
                    simu_name = 'simu' + simu_id
                    newgroup = db.createGroup(db.root, simu_name)
                    f.copyNode(f.root, newgroup, recursive=True)
                    f.copyNodeAttrs(f.root, newgroup)


def get_all_simus_id(dbpath):
    """Return simulations id (time, uuid) from the DB"""
    ids = []
    with tables.openFile(dbpath) as db:
        for g in db.walkGroups():
            try:
                ids.append((g._v_attrs['time'], g._v_attrs['uuid']))
            except KeyError:
                pass
    return ids
