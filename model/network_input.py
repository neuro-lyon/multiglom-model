#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import Equations
from brian.units import *
from brian.stdunits import *
from parameters import Input as InputParameters

class NetworkInput:
    """Input for the network of mitral and granule cells."""

    # Parameters are global so Brian can put them into the equations.
    # You can change their values in parameters.py
    global tau_Ein, g_Ein0, sigma_Ein
    psin = InputParameters()
    tau_Ein   = psin.tau_Ein
    g_Ein0    = psin.g_Ein0
    sigma_Ein = psin.sigma_Ein

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
