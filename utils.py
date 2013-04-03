"""
Utilities not related to neuronal network simulation.

"""

import importlib
import model

def path_to_modline(filepath):
    """Returns a string representation for module import of a file."""
    # Import the set of parameters
    pathsplit = filepath.split('/')
    modname = pathsplit[-1][:-3]
    return '.'.join(pathsplit[:-1] + [modname])

def set_model_ps(filepath, dicname='PARAMETERS'):
    """Set the model parameters given the parameter set in `filepath`."""
    psmod = importlib.import_module(path_to_modline(filepath))
    model.PARAMETERS = getattr(psmod, dicname)
