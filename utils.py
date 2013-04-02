"""
Utilities not related to neuronal network simulation.

"""

from os import getcwd, chdir

def get_module_var(filepath, varname):
    """Returns the variable `varname` inside the file `filepath`."""
    fromdir = getcwd()
    # Import the set of parameters
    pathsplit = filepath.split('/')
    directory = '/'.join(pathsplit[:-1])
    modname = pathsplit[-1][:-3]
    chdir(directory)
    mod = __import__(modname)
    var = getattr(mod, varname)
    chdir(fromdir)
    return var
