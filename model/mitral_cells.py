#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import Equations, NeuronGroup
from brian.units import *
from brian.stdunits import *
from parameters import Mitral as MitralParameters

class MitralCells:
    """Population of mitral cells."""

    # Parameters are global so Brian can put them into the equations.
    # You can change their values in parameters.py
    global N, V_t, V_r, t_refract
    global C_m, g_L, E_L
    psmt = MitralParameters()
    V_t = psmt.V_t
    V_r = psmt.V_r
    t_refract = psmt.t_refract
    C_m  = psmt.C_m
    g_L = psmt.g_L
    E_L = psmt.E_L

    def __init__(self):
        """Create an empty population of mitral cells."""
        self.eqs_model = Equations()
        self.pop = None

    def add_eqs(self, supp_eqs=None):
        """Add the standard equation of the LIF model.

        Supplementary equations are added through `supp_eqs` this way:
        supp_eqs = {'var': ['-I_input', '+I_electrode'],
                    'eqs': [Equations('I_input = g_Ein*(V - 0*mvolt) : amp*meter**-2'),
                            Equations('I_electrode : amp*meter**-2')]}
        """
        t = ''
        # Add the extra variable into the LIF equation
        if supp_eqs:
            t = ' '.join([vname for vname in supp_eqs['var']])
            # Add the extra equations
            for eq in supp_eqs['eqs']:
                self.eqs_model += eq
        std_lif_eq = Equations('dV/dt = (-g_L*(V - E_L) '+ t +')/C_m : volt')
        self.eqs_model += std_lif_eq

    def make_pop(self, N):
        """Makes the NeuronGroup for a population of `N` mitral cells.

        This must be called after all the other construction methods.
        """
        self.pop = NeuronGroup(N, model=self.eqs_model, order=1,
                               threshold=V_t, reset=V_r, refractory=t_refract)
