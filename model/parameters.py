#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from sys import modules
from os import chdir


# Set the parameter file here (path relative to python call)
PARAMETER_FILE = 'paramsets/std.py'


class Empty():
    pass


def set_params(filepath=PARAMETER_FILE, dic='parameters'):
    # Import the set of parameters
    pathsplit = filepath.split('/')
    directory = '/'.join(pathsplit[:-1])
    modname = pathsplit[-1][:-3]
    chdir(directory)
    mod = __import__(modname)
    paramset = getattr(mod, dic)
    # Build the parameter classes and attributes accordingly
    for classname in paramset:
        setattr(modules[__name__], classname, Empty)
        for parname in paramset[classname]:
            curclass = getattr(modules[__name__], classname)
            setattr(curclass, parname, paramset[classname][parname])


set_params()
