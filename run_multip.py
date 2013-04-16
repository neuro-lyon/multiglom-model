# -*- coding:utf-8 -*-
"""
Script for running multiple simulation with different parameter set.

"""
from multiprocessing import Pool
from os import listdir, path
from sys import argv


def new_simu(psfile):
    # Creating simulation run arguments
    from arg_parser import APARSER
    args = APARSER.parse_args()
    args.no_plot = True


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
    n_processes = int(argv[1])
    dir_paramsets = argv[2]

    # Get all the parameter set files
    psfiles = get_pset_files(dir_paramsets)

    # Run the simulation with each parameter set
    pool = Pool(processes=n_processes)
    pool.map(run_simu, psfiles)
