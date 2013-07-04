# -*- coding:utf-8 -*-

from brian import *

from model import PARAMETERS as ps

PSIN = ps['Input']
tau_Ein   = PSIN['tau_Ein']
g_Ein0    = PSIN['g_Ein0']
sigma_Ein = PSIN['sigma_Ein']

PSINOSC = ps['InputOscillation']
f    = PSINOSC['f']
C    = PSINOSC['C']


class Glomerule:
    """Ensemble of glomeruli, driving input to mitral cells."""

    def __init__(self):
        """Create an empty glomerule."""
        self.eqs_model = Equations()
        self.pop = None
    
    def add_eqs(self, oscillating=False):
        if oscillating:
            std_eqs = Equations('''
                dn/dt = -(n + tau_Ein**.5*xi)/tau_Ein : 1
                sigma = sigma_Ein*sqrt(1 + C*cos(2*pi*f*t)) : siemens*meter**-2
                m0 = g_Ein0*(1 + C*cos(2*pi*f*t)) : siemens*meter**-2
                dm/dt = (m0 - m)/tau_Ein : siemens*meter**-2
                g = m + sigma*n: siemens*meter**-2
                ''')
        else:
            std_eqs = Equations('''
                dg/dt = (g_Ein0 - g + sigma_Ein*tau_Ein**.5*xi)/tau_Ein : siemens*meter**-2
                ''')
        self.eqs_model += std_eqs

    def make_pop(self, N):
        self.pop = NeuronGroup(N, model=self.eqs_model)
