#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from sys import modules
from os import chdir, getcwd


# Set the parameter file here (path relative to python call)
PARAMETER_FILE = 'paramsets/std_beta.py'


class MakeupClass():
    def __init__(self):
        pass

    def __call__(self):
        return self

    def add_attr(self, name, value):
        setattr(self, name, value)


def set_params(filepath=PARAMETER_FILE, dic='parameters'):
    fromdir = getcwd()
    # Import the set of parameters
    pathsplit = filepath.split('/')
    directory = '/'.join(pathsplit[:-1])
    modname = pathsplit[-1][:-3]
    chdir(directory)
    mod = __import__(modname)
    chdir(fromdir)
    paramset = getattr(mod, dic)
    # Build the parameter classes and attributes accordingly
    for classname in paramset:
        curclass = MakeupClass()
        setattr(modules[__name__], classname, curclass)
        for parname in paramset[classname]:
            curclass.add_attr(parname, paramset[classname][parname])


set_params()
