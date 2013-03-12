#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import Equations, NeuronGroup
from brian.units import *
from brian.stdunits import *
from parameters import Granule as GranuleParameters

class GranuleCells:
    """Population of granule cells."""

    def __init__(self):
        """Create an empty population of granule cells."""
        self.eqs_model = Equations()
        self.pop = None

        psgr = GranuleParameters()
        self.C_m  = psgr.C_m
        self.g_L  = psgr.g_L
        self.E_L  = psgr.E_L
        self.g_SD = psgr.g_SD
        self.g_DS = psgr.g_DS

    def add_eqs(self, supp_eqs=None):
        """Add the standard equations of granule model.

        Supplementary equations are added to `V_D`through `supp_eqs` this way:
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
        std_eqs  = Equations('dV_D/dt = (-g_L*(V_D - E_L) + g_DS*(V_S - V_D) '+ t +')/C_m : mvolt')
        std_eqs += Equations('dV_S/dt = (-g_L*(V_S - E_L) + g_SD*(V_D - V_S))/C_m : mvolt')
        self.eqs_model += std_eqs

    def make_pop(self, N):
        """Makes the NeuronGroup for a population of `N` granule cells.

        This must be called after all the other construction methods.
        """
        self.pop = NeuronGroup(N, model=self.eqs_model, order=1, implicit=True)

