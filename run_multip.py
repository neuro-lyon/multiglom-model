# -*- coding:utf-8 -*-
"""
Script for running multiple simulation with different parameter set.

"""
from subprocess import check_output
from multiprocessing import Pool
from os import listdir, mkdir, path
from sys import argv


def run_simu((psfile, dir_runs)):
    """Runs the main script with specified parameter file."""
    ndir = psfile[:-3]
    mkdir(ndir)
    simu_output = check_output(["python2.7", "multiglom_network.py", psfile,
                                "--full-ps"])
    with open(path.join(ndir, "output.txt"), 'w') as f:
        f.write(simu_output)


if __name__ == '__main__':
    n_processes = int(argv[1])
    dir_runs = argv[2]

    psfiles = listdir(dir_runs)
    psfile_filter = lambda fname: fname[-3:] == '.py' and fname != '__init__.py'
    psfiles = filter(psfile_filter, psfiles)
    psfiles = [(path.join(dir_runs, fname), dir_runs) for fname in psfiles]

    pool = Pool(processes=n_processes)
    pool.map(run_simu, psfiles)
