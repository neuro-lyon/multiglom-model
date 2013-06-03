# -*- coding:utf-8 -*-

import tables
import os
from utils import listdir_filter


def init_data_h5(filename):
    """Initialize a data HDF5 file"""
    filename = filename_to_h5(filename)
    if not file_exists(filename):
        f = tables.openFile(filename, 'w')
        f.close()


def write_simu_data(filename, info, parameters, results):
    """Create a HDF5 file for new simulation"""
    with tables.openFile(filename, 'a') as f:
        # Put info into the HDF5 root
        for attr in info:
            setattr(f.root._v_attrs, attr, info[attr])

        # Put the parameter set into the file
        ps = f.createGroup('/', 'paramset', title="Parameters")
        for attr in parameters['set']:
            setattr(ps._v_attrs, attr, parameters['set'][attr])
        ps_arrays = f.createGroup(ps, 'arrays', title="Parameter arrays")
        for ps_array in parameters['arrays']:
            f.createArray(ps_arrays, ps_array, parameters['arrays'][ps_array][0],
                          title=parameters['arrays'][ps_array][1])

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


def collect_h5_to_db(dirpath, dbpath, output=True):
    """Collects all HDF5 in dirpath and put them into the big HDF5 DB"""
    h5files = listdir_filter(dirpath, lambda fname: fname[-3:] == '.h5')
    try:
        h5files.remove(dbpath)  # In case db is in the same directory
    except ValueError:
        pass  # Nothing to do if the db path is not in the same directory

    n_files = len(h5files)
    old_simus_id = get_all_simus_id(dbpath)

    # Open the db to put the new simus in
    with tables.openFile(dbpath, 'a') as db:
        for num_file, h5file in enumerate(h5files):
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

            # Output progress
            if output:
                print 'File '+str(num_file + 1)+'/'+str(n_files)+' done.'


def get_all_simus_id(dbpath):
    """Return simulations id (time, uuid) from the DB"""
    ids = []
    with tables.openFile(dbpath) as db:
        for g in get_first_level_groups(db.root):
            ids.append((g._v_attrs['time'], g._v_attrs['uuid']))
    return ids


def get_first_level_groups(root_group):
    """Return an non-ordered list of first-level groups of root_group"""
    groups = []
    for group_name in root_group._v_groups:
        groups.append(getattr(root_group, group_name))
    return groups


def get_all_attrs(db, attrs_path):
    """Get the attributes for each group of db (except '/')

    Parameters
    ----------
    db : tables.file.File
    attrs_path : 2-D list
    """
    groups_attrs = []
    for group in get_first_level_groups(db.root):
        attrs = []
        for apath in attrs_path:
            attrs.append(get_group_attr(group, apath))
        groups_attrs.append(attrs)
    return groups_attrs


def get_group_attr(group, attr_path):
    """Return the attribute specified by a path in a group"""
    res = group
    for attr in attr_path:
        if type(res) == type({}):
            res = res.get(attr)
        else:
            try:
                res = getattr(res, attr)
            except tables.NoSuchNodeError, e:
                print e
    return res
