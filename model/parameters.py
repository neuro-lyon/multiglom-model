#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian.units import *
from brian.stdunits import *

class Common:
    """General parameters."""
    def __init__(self):
        self.simu_dt     = 0.05*msecond
        self.simu_length = 2000*msecond
        self.N_subpop    = 1
        self.N_mitral    = 100

class Glomerule:
    """Parameters for the glomeruli."""
    def __init__(self):
        self.tau = 3*msecond
        self.f   = 2*Hz
        self.A   = .1*second**-.5*siemens*meter**-2
        self.B   = 10*second**-1*siemens*meter**-2
        self.C   = 1

class Mitral:
    """Parameters for the mitral cell population."""
    def __init__(self):
        """Define parameter values here."""
        self.C_m       = 0.08*farad*meter**-2
        self.g_L       = 0.87*siemens*meter**-2
        self.E_L       = -64.5*mvolt
        self.V_r       = -74*mvolt
        self.V_t       = -62*mvolt
        self.t_refract = 0.2*msecond

class Granule:
    """Parameters for the granule cell population."""
    def __init__(self):
        self.C_m  = 0.01*farad*meter**-2
        self.g_L  = 0.83*siemens*meter**-2
        self.E_L  = -70*mvolt
        self.g_SD = 1*siemens*meter**-2
        self.g_DS = 300*siemens*meter**-2

class Input:
    """Input parameters."""
    def __init__(self):
        self.tau_Ein   = 3*msecond
        self.g_Ein0    = 0.6*siemens*meter**-2
        sigma_Ein_low  = .05*siemens*meter**-2*second**(-1./2)
        sigma_Ein_high = 5.*sigma_Ein_low
        self.sigma_Ein = sigma_Ein_low

class Synapse:
    """Synapse parameters."""
    def __init__(self):
        self.V_E      = 0*mvolt
        self.V_act_E  = 0*mvolt
        g_E_low       = 0.7*siemens*meter**-2
        g_E_high      = 3.5*siemens*meter**-2
        self.g_E      = g_E_high
        self.sigma_E  = 0.01*mvolt
        self.alpha_E  = 10*msecond**-1
        self.beta_E   = 1./3*msecond**-1

        self.V_I      = -70*mvolt
        self.V_act_I  = -66.4*mvolt
        self.g_I      = 10*siemens*meter**-2
        self.sigma_I  = 0.4*mvolt
        self.alpha_I  = 5*msecond**-1
        self.beta_I   = 1./10*msecond**-1
