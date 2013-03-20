#!/usr/bin/env python2
# -*- coding:utf-8 -*-

from brian import *
from parameters import Glomerule as GlomeruleParameters
from parameters import Input as InputParameters

psin = InputParameters()
tau_Ein   = psin.tau_Ein
g_Ein0    = psin.g_Ein0
sigma_Ein = psin.sigma_Ein

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
    
    def add_eqs(self, oscillation=True):
        if oscillation:
            std_eqs = Equations('''
                dn/dt = -(n + xi)/tau : second**-.5
                sigma = A*tau**.5*sqrt(1 + C*cos(2*pi*f*t)) : siemens*meter**-2
                m0 = B*tau*(1 + C*cos(2*pi*f*t)) : siemens*meter**-2
                dm/dt = (m0 - m)/tau : siemens*meter**-2
                g = m + sigma*n*tau**.5 : siemens*meter**-2
                ''')
        else:
            std_eqs = Equations('''
                dg/dt = (g_Ein0 - g)/tau_Ein + sigma_Ein * xi : siemens*meter**-2
                ''')
        self.eqs_model += std_eqs

    def make_pop(self, N):
        self.pop = NeuronGroup(N, model=self.eqs_model)

if __name__ == '__main__':
    clear(erase=True, all=True)
    defaultclock.reinit()

    g = Glomerule()
    nglom = 3
    g.add_eqs()
    g.make_pop(nglom)

    m = StateMonitor(g.pop, 'g', record=True)
    netw = Network(g.pop, m)
    netw.run(500*msecond)

    for i in xrange(nglom):
        plot(m.times/msecond, m[i]/(siemens*meter**-2))
    show()
