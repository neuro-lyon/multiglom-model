#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import *
from parameters import Glomerule as GlomeruleParameters

psgm = GlomeruleParameters()
tau  = psgm.tau
f    = psgm.f
A    = psgm.A
B    = psgm.B
C    = psgm.C

class Glomerule:
    """Ensemble of glomeruli, driving input to mitral cells."""

    def __init__(self):
        """Create an empty glomerule."""
        self.eqs_model = Equations()
        self.pop = None
    
    def add_eqs(self):
        std_eqs = Equations('''
            dn/dt = -(n + xi)/tau : second**-.5
            sigma = A*tau**.5*sqrt(1 + C*cos(2*pi*f*t)) : 1
            m0 = B*tau*(1 + C*cos(2*pi*f*t)) : 1
            dm/dt = (m0 - m)/tau : 1
            g = m + sigma*tau**.5*n : 1
        ''')
        self.eqs_model += std_eqs

    def make_pop(self, N):
        self.pop = NeuronGroup(N, model=self.eqs_model)
