# -*- coding:utf-8 -*-

from brian import Equations
from brian.units import *
from brian.stdunits import *

from model import PARAMETERS as ps

PSIN = ps['Input']
tau_Ein   = PSIN['tau_Ein']
g_Ein0    = PSIN['g_Ein0']
sigma_Ein = PSIN['sigma_Ein']

class NetworkInput:
    """Input for the network of mitral and granule cells."""

    def __init__(self):
        """Create an empty model of equations defining the network input"""
        self.eqs_model = Equations()

    def get_eqs_model(self):
        """Get the model of equations"""
        return self.eqs_model

    def set_eqs_model(self, eqs=None):
        """Sets the model of equations.
        If no model if specified, uses the standard model.
        """
        if eqs:
            self.eqs_model = eqs
        else:
            self.eqs_model = Equations('''
            dg_Ein/dt = (g_Ein0 - g_Ein)/tau_Ein + sigma_Ein * xi : siemens*meter**-2
            I_input = g_Ein*(V - 0*mvolt) : amp*meter**-2
            ''')
