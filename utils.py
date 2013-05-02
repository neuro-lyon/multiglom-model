# -*- coding:utf-8 -*-
"""
Utilities not related to neuronal network simulation.

"""
from numpy import linspace
from os import path, listdir
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
    params = {('Input', 'g_Ein0'): {'range': linspace(start=0, stop=1, num=4),
                                    'unit': siemens*meter**-2},
              ('Common', 'inter_conn_rate', '*', '*'): {'range': linspace(start=0.2, stop=0.8, num=2),
                                                        'unit': 1}
    }
    """
    template = get_template(template_file)

    # Map params dict key to indices
    params_map = []
    for key in params:
        params_map.append(key)

    # Build range list
    range_list = []
    for key, value in zip(params.keys(), params.values()):
        range_list.append(value['range'])

    # Iterate on the cartesian product of the ranges to create new sets
    index=0
    for comb in itertools.product(*range_list):
        index+=1
        fname = 'parset'+str(0)*(5-len(str(index)))+str(index)+"__" # assume less than 10000 parameters sets

        # Set the new value for each variable of the combination
        for ind_var in xrange(len(comb)):
            var_key = params_map[ind_var]
            var_value = comb[ind_var]
            var_unit = params[var_key]['unit']

            change_dict_key(template, var_key, var_value*var_unit)

            fname += '__'.join([_.replace('*', '') for _ in var_key])
            fname += '_' + str(comb[ind_var]).replace('.', '_') + '_'
        fname += '.py'

        # Write the new parameter set to a file
        with open(path.join(output_dir, fname), 'w') as f:
            f.writelines(["from brian.stdunits import *\n",
                          "from brian.units import *\n"])
            f.write('PARAMETERS = ' + str(template) + '\n')
            

    # Put a __init__.py to make the modules importable
    f = open(path.join(output_dir, '__init__.py'), 'w')
    f.close()


def get_template(template_file, varname="PARAMETERS"):
    """Return the parameter variable from a template_file"""
    module = import_module(path_to_modline(template_file))
    return getattr(module, varname)


def change_dict_key(dic, path, new_value, anykey='*'):
    """Change0 *in place* the dictionnary key value to new_value.

    If key is '*', then change all dictionnary key values to new_value.
    """
    # Stop case
    if len(path) == 1:
        key = path[0]
        if key == anykey:
            for k in dic:
                dic[k] = new_value
        elif dic.has_key(key):
            dic[key] = new_value

    # Recursion
    else:
        sub = path[0]
        if sub == anykey:
            for k in dic:
                change_dict_key(dic[k], path[1:], new_value, anykey)
        else:
            change_dict_key(dic[sub], path[1:], new_value, anykey)

    return dic


def pairs(n, no_ident=True):
    """All possible pairs (i, j) with i and j < n"""
    list_pairs = []
    for i in xrange(n):
        for j in xrange(i, n):
            if not (no_ident and i == j):
                list_pairs.append((i, j))
    return list_pairs


def listdir_filter(dir, file_filter):
    """List all files in dir and apply a filter to select them"""
    files = listdir(dir)
    filtered_files = filter(file_filter, files)
    return [path.join(dir, fname) for fname in filtered_files]


def to1d(x, y, lenx):
    """1-D coordinate from 2-D coordinate"""
    return x + y*lenx


from brian import *
if __name__ == '__main__':
    d = {('Common', 'inter_conn_strength', '*', '*'):
            {'range': linspace(start=0., stop=1., num=10),
             'unit': 1},
         ('Common', 'inter_conn_rate', '*', '*'):
            {'range': linspace(start=0., stop=1, num=10),
             'unit': 1}
    }
    gen_parameters('paramsets/std_beta.py', d, 'runs/interco10x10/beta')
