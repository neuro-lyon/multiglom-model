# -*- coding:utf-8 -*-
"""
Script for running multiple simulation with different parameter set.

"""
from multiprocessing import Pool
from os import listdir, path
from datetime import datetime
from uuid import uuid4
import git

from arg_parsers import SIM_PARSER, MULTISIM_PARSER
import multiglom_network
from data_collection.h5manager import init_data_h5, write_simu_data


def new_simu(psfile):
    """Run a simulation using the specified parameter set file"""
    filedir = path.dirname(psfile) +  '/'
    # Register system state
    info = get_sys_state()

    # Creating simulation run arguments
    args = SIM_PARSER.parse_args([psfile, '--no-plot'])

    # Run simulation and get results
    paramset, results = multiglom_network.main(args)

    # Write the simulation to HDF5
    nfilename = filedir + info['time'].replace(':', '-')
    nfilename += '_' + info['uuid'] + '.h5'
    init_data_h5(nfilename)
    write_simu_data(nfilename, info, paramset, results)


def get_sys_state():
    """Get system state"""
    sys_state = {}
    # Time and UUID
    sys_state['time'] = str(datetime.now()).replace(' ', '_')
    sys_state['uuid'] = str(uuid4())

    # Git revision
    repo = git.Repo('.')
    assert not repo.is_dirty(), "Your git repository has uncommited changes."
    sys_state['git_rev'] = str(repo.commit(repo.head))

    return sys_state


def get_pset_files(dir_paramsets):
    """Return a list of parameter set files"""
    psfiles = listdir(dir_paramsets)

    def pset_filter(fname):
        return fname[-3:] == '.py' and fname != '__init__.py'
    psfiles = filter(pset_filter, psfiles)

    psfiles = [path.join(dir_paramsets, fname) for fname in psfiles]
    return psfiles


if __name__ == '__main__':
    # Get run arguments
    args = MULTISIM_PARSER.parse_args()

    # Get all the parameter set files
    psfiles = get_pset_files(args.pset_dir)

    # Run the simulation with each parameter set
    pool = Pool(processes=args.nproc)
    pool.map(new_simu, psfiles)
