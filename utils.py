"""
Utilities not related to neuronal network simulation.

"""

def path_to_modline(filepath):
    """Returns a string representation for module import of a file."""
    # Import the set of parameters
    pathsplit = filepath.split('/')
    modname = pathsplit[-1][:-3]
    return '.'.join(pathsplit[:-1] + [modname])
