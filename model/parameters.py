#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian.units import *
from brian.stdunits import *

class Mitral:
    """Parameters for the mitral cell population."""
    def __init__(self):
        """Define parameter values here."""
        self.C_m       = 0.08*farad*meter**-2
        self.g_L       = 0.87*siemens*meter**-2
        self.E_L       = -64.5*mvolt
        # self.V_r       = -69.8*mvolt # from tests to match F. model
        # self.V_t       = -58*mvolt # from tests to match F. model
        self.V_r       = -74*mvolt
        self.V_t       = -62*mvolt
        self.t_refract = 0.20000001*msecond # pour éviter le bug effet de bord

class Granule:
    """Parameters for the granule cell population."""
    def __init__(self):
        self.C_m  = 0.01*farad*meter**-2
        self.g_L  = 0.83*siemens*meter**-2
        self.E_L  = -70*mvolt
        self.g_SD = 1*siemens*meter**-2
        self.g_DS = 300*siemens*meter**-2

class Input:
    def __init__(self):
        self.tau_Ein   = 3*msecond
        self.g_Ein0    = 0.5*siemens*meter**-2
        self.sigma_Ein = 2*siemens*meter**-2*second**(-1./2)

class Synapse:
    def __init__(self):
        self.V_E     = 0*mvolt
        self.V_act_E = 0*mvolt
        g_E_low      = 0.7*siemens*meter**-2
        g_E_high     = 3.5*siemens*meter**-2
        self.g_E     = g_E_high
        self.sigma_E = 0.01*mvolt
        self.alpha_E = 10*msecond**-1
        self.beta_E  = 1./3*msecond**-1

        self.V_I      = -70*mvolt
        self.V_act_I  = -66.4*mvolt
        self.g_I      = 10*siemens*meter**-2
        self.sigma_I  = 0.4*mvolt
        self.alpha_I  = 5*msecond**-1
        self.beta_I   = 1./10*msecond**-1
