# -*- coding:utf-8 -*-
"""
Utilities not related to neuronal network simulation.

"""
from numpy import arange
from os import path
from copy import deepcopy
from importlib import import_module
from itertools import product

def path_to_modline(filepath):
    """Returns a string representation for module import of a file."""
    # Import the set of parameters
    pathsplit = filepath.split('/')
    modname = pathsplit[-1][:-3]
    return '.'.join(pathsplit[:-1] + [modname])


def print_dict(dictio, level=0):
    """Pretty print a dictionnary."""
    for k in dictio:
        if type(dictio[k]) == type({}):
            print ' '*2*level + str(k)
            print_dict(dictio[k], level + 1)
        else:
            print ' '*2*level + str(k) + ': ' + str(dictio[k])


def gen_parameters(template_file, params, output_dir):
    """Generate parameter sets from template.

    template: filename with a dictionnary of parameters
    params: dictionnary of parameters with start, stop, step and unit
    output_dir: directory to put the created parameter sets in

    Example of params:
    params = {'Input': {'g_Ein0': {'start': 0.,
                                  'stop' : 1.,
                                  'step' : 0.2,
                                  'unit' : siemens*meter**-2}}}
    """
    original_template = get_template(template_file)
    for category in params:
        for var in params[category]:
            # Get the info for the parameter to change
            var_details = params[category][var]
            var_range = arange(var_details['start'], var_details['stop'], var_details['step'])
            var_units = var_details['unit']

            # Create a new file for each variable value
            for value in var_range:
                new_parameters = deepcopy(original_template)
                new_parameters[category][var] = value*var_units
                units = str(var_units).replace(' ', '_')
                fname = '__'.join([category, var, str(value), units])
                fname = fname.replace('.', '_')  # So we can import it easily as a module
                fname += '.py'
                with open(path.join(output_dir, fname), 'w') as f:
                    f.write('PARAMETERS = ' + str(new_parameters))

    # Put a __init__.py to make the modules importable
    f = open(path.join(output_dir, '__init__.py'), 'w')
    f.close()


def get_template(template_file, varname="PARAMETERS"):
    """Return the parameter variable from a template_file"""
    module = import_module(path_to_modline(template_file))
    return getattr(module, varname)


from brian import *
if __name__ == '__main__':
    d = {"Input":
            {'g_Ein0': {'start': 0., 'stop': 1, 'step': 0.5, 'unit': siemens*meter**-2},
            'tau_Ein': {'start': 0., 'stop': 3, 'step': 0.5, 'unit': second}}
        }
    gen_parameters('paramsets/std_beta.py', d, '/tmp/ps2')
