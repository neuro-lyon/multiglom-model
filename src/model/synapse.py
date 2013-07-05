# -*- coding:utf-8 -*-

from brian import Equations
from brian.units import *
from brian.stdunits import *

from model import PARAMETERS as ps

# If the model parameters have been initialized
if ps != None:
    N_subpop = ps['Common']['N_subpop']

    PSSYN = ps['Synapse']
    V_E      = PSSYN['V_E']
    V_act_E  = PSSYN['V_act_E']
    g_E      = PSSYN['g_E']
    sigma_E  = PSSYN['sigma_E']
    alpha_E  = PSSYN['alpha_E']
    beta_E   = PSSYN['beta_E']

    V_I      = PSSYN['V_I']
    V_act_I  = PSSYN['V_act_I']
    g_I      = PSSYN['g_I']
    sigma_I  = PSSYN['sigma_I']
    alpha_I  = PSSYN['alpha_I']
    beta_I   = PSSYN['beta_I']


class Synapse:
    """A Synapse is a connection between a mitral cell and a granule cell.

    Synapses are modeled with a set of ordinary differential equations (ODE).
    There are two kind of synapses: inhibitory and excitatory. Each kind of
    synapse has its set of ODE, but they share common features.

    Attributes
    ----------
    is_exc : bool
        True if the synapse is excitatory, False if inhibitory
    is_inhib : bool
        True if the synapse is inhibitory, False if excitatory
    eqs_model : brian.Equations
        Set of equations that define the model

    """

    def __init__(self, synapse_type):
        """Initialize a Synapse.

        At this point, the equation model is empty. The user must define the
        equation model after the initialization by using the
        :meth:`set_eqs_model`.
        """
        if synapse_type[:3] == 'exc':
            self.is_exc = True
        elif synapse_type[:3] == 'inh':
            self.is_exc = False
        self.is_inhib = not self.is_exc
        self.eqs_model = Equations()

    def get_eqs_model(self):
        """Get the model of equations of this synapse"""
        return self.eqs_model

    def set_eqs_model(self, eqs=None):
        """Sets the model of equations.

        A standard model is used if no model is specified with :attr:`eqs`.
        If no model if specified, uses the standard model (according to the
        synapse type).

        Parameters
        ----------
        eqs : brian.Equations
            Set of equations that define the model

        Notes
        -----
        This method must be called after creating a :class:`Synapse`, otherwise
        the equation model is empty.
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
                    s_syn_self : 1
                ''')
            else:
                self.eqs_model = Equations('''
                    I_syn = g_I * s_syn * (V - V_I) : amp*meter**-2
                    ds/dt = alpha_E * (1 - s) * T - beta_E * s : 1
                    T : 1
                    s_syn : 1
                ''')
