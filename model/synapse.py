#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import Equations
from brian.units import *
from brian.stdunits import *
from parameters import Synapse as SynapseParameters

class Synapse:
    """Synapse, from mitral cells to granule cells."""

    # Parameters are global so Brian can put them into the equations.
    # You can change their values in parameters.py
    global V_E, V_act_E, g_E, sigma_E, alpha_E, beta_E, V_E
    global V_I, V_act_I, g_I, sigma_I, alpha_I, beta_I, V_I
    pssyn = SynapseParameters()
    V_E      = pssyn.V_E 
    V_act_E  = pssyn.V_act_E
    g_E      = pssyn.g_E
    sigma_E  = pssyn.sigma_E
    alpha_E  = pssyn.alpha_E
    beta_E   = pssyn.beta_E

    V_I      = pssyn.V_I
    V_act_I  = pssyn.V_act_I
    g_I      = pssyn.g_I
    sigma_I  = pssyn.sigma_I
    alpha_I  = pssyn.alpha_I
    beta_I   = pssyn.beta_I

    def __init__(self, synapse_type):
        if synapse_type[:3] == 'exc':
            self.is_exc = True
        elif synapse_type[:3] == 'inh':
            self.is_exc = False
        self.is_inhib = not self.is_exc
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
            if self.is_inhib:
                self.eqs_model = Equations('''
                    I_syn = g_E * s_syn * (V_D - V_E) : amp*meter**-2
                    ds/dt = alpha_I * (1 - s) * T - beta_I * s : 1
                    T = 1/(1 + exp(-1*(V_D - V_act_I)/sigma_I)) : 1
                    s_syn : 1
                ''')
            else:
                self.eqs_model = Equations('''
                    I_syn = g_I * s_syn * (V - V_I) : amp*meter**-2
                    ds/dt = alpha_E * (1 - s) * T - beta_E * s : 1
                    T : 1
                    s_syn : 1
                ''')
