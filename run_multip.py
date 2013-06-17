# -*- coding:utf-8 -*-
"""
Script for running multiple simulation with different parameter set.

"""
from multiprocessing import Pool
from os import path
from datetime import datetime
from uuid import uuid4
from math import ceil
import git

from arg_parsers import SIM_PARSER, MULTISIM_PARSER
import multiglom_network
from h5manager import init_data_h5, write_simu_data
from utils import listdir_filter


def new_simu((psfile, num_file, tot_file, check_repo_is_dirty)):
    """Run a simulation using the specified parameter set file"""
    filedir = path.dirname(psfile) +  '/'
    # Register system state
    info = get_sys_state(check_repo_is_dirty)

    # Creating simulation run arguments
    args = SIM_PARSER.parse_args([psfile, '--no-plot', '--no-brian-output',
                                  '--no-summary'])

    # Run simulation and get results
    paramset, results = multiglom_network.main(args)

    # Write the simulation to HDF5
    nfilename = filedir + info['time'].replace(':', '-')
    nfilename += '_' + info['uuid'] + '.h5'
    init_data_h5(nfilename)
    write_simu_data(nfilename, info, paramset, results)


def get_sys_state(check_repo_is_dirty=True):
    """Get system state"""
    sys_state = {}
    # Time and UUID
    sys_state['time'] = str(datetime.now()).replace(' ', '_')
    sys_state['uuid'] = str(uuid4())

    # Git revision
    repo = git.Repo('.')

    is_repo_dirty = repo.is_dirty()
    if check_repo_is_dirty:
        assert not is_repo_dirty, \
            "Your git repository has uncommited changes, use --testing for quick testing without asserting that the repo is clean."
    sys_state['git_repo_is_dirty'] = is_repo_dirty

    sys_state['git_rev'] = str(repo.commit(repo.head))

    return sys_state


def get_pset_files(dir_paramsets):
    """Return a list of parameter set files"""
    def pset_filter(fname):
        return fname[-3:] == '.py' and fname != '__init__.py'

    return listdir_filter(dir_paramsets, pset_filter)


if __name__ == '__main__':
    # Get run arguments
    args = MULTISIM_PARSER.parse_args()

    # Get all the parameter set files
    psfiles = get_pset_files(args.pset_dir)
    check_repo_is_dirty = not args.testing

    # Run the simulation with each parameter set
    simu_args = []
    for ind, psfile in enumerate(psfiles):
        simu_args.append((psfile, ind, len(psfiles), check_repo_is_dirty))

    print "Start simul"
    npool = int(ceil(1.*len(simu_args)/args.nproc))
    for i in range(npool):
        pool = Pool(processes=args.nproc)
        try:
            pool.map(new_simu, simu_args[i*args.nproc:(i+1)*args.nproc])
        except Exception, e:
            print "WARNING:", e
            print "SIMULATION %s TO %s MIGHT BE INVALID" % (i*args.nproc,
                                                           (i+1)*args.nproc)
        pool.close()
    print "End simul"
