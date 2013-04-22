# -*- coding:utf-8 -*-
"""
Utilities not related to neuronal network simulation.

"""
from numpy import arange
from os import path
from copy import deepcopy
from importlib import import_module
import itertools

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
    """Generate parameter sets from a template.

    template: filename with a dictionnary of parameters in it
    params: dictionnary of parameters with start, stop, step and unit
    output_dir: directory to put the created parameter sets in

    Example of params:
    params = {"Input":
                {'g_Ein0': {'start': 0.,
                            'stop': 1,
                            'step': 0.5,
                            'unit': siemens*meter**-2},
                 'tau_Ein': {'start': 0.,
                             'stop': 3,
                             'step': 0.5,
                             'unit': second}}}
    """
    template = get_template(template_file)
    var_list = []
    range_list = []

    # Build the variable value ranges
    for category in params:
        for var in params[category]:
            # Get the info for the parameter to change
            var_details = params[category][var]
            var_range = arange(var_details['start'], var_details['stop'], var_details['step'])
            var_units = var_details['unit']

            var_list.append({'name': var, 'units': var_units, 'cat': category})
            range_list.append(var_range)

    # Iterate on the cartesian product of the ranges to create new sets
    for comb in itertools.product(*range_list):
        fname = ''
        for ind_var in xrange(len(var_list)):
            var_category = var_list[ind_var]['cat']
            var_name = var_list[ind_var]['name']
            var_units = var_list[ind_var]['units']
            template[var_category][var_name] = comb[ind_var]*var_units
            var_value = str(comb[ind_var]).replace('.', '_')
            fname += '__'.join([var_category, var_name, var_value])
        fname += '.py'
        with open(path.join(output_dir, fname), 'w') as f:
            f.writelines(["from brian.stdunits import *",
                          "from brian.units import *", "\n"])
            f.write('PARAMETERS = ' + str(template) + '\n')

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
             'tau_Ein': {'start': 0., 'stop': 3, 'step': 1, 'unit': second}
            }
        }
    gen_parameters('paramsets/std_beta.py', d, '/tmp/ps')
