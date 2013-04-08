# -*- coding:utf-8 -*-

from brian import *

from model import PARAMETERS as ps

PSIN = ps['Input']
tau_Ein   = PSIN['tau_Ein']
g_Ein0    = PSIN['g_Ein0']
sigma_Ein = PSIN['sigma_Ein']

PSGM = ps['Glomerule'] 
tau  = PSGM['tau']
f    = PSGM['f']
A    = PSGM['A']
B    = PSGM['B']
C    = PSGM['C']

class Glomerule:
    """Ensemble of glomeruli, driving input to mitral cells."""

    def __init__(self):
        """Create an empty glomerule."""
        self.eqs_model = Equations()
        self.pop = None
    
    def add_eqs(self, oscillating=True):
        if oscillating:
            std_eqs = Equations('''
                dn/dt = -(n + xi)/tau : second**-.5
                sigma = A*tau**.5*sqrt(1 + C*cos(2*pi*f*t)) : siemens*meter**-2
                m0 = B*tau*(1 + C*cos(2*pi*f*t)) : siemens*meter**-2
                dm/dt = (m0 - m)/tau : siemens*meter**-2
                g = m + sigma*n*tau**.5 : siemens*meter**-2
                ''')
        else:
            std_eqs = Equations('''
                dg/dt = (g_Ein0 - g)/tau_Ein + sigma_Ein*xi : siemens*meter**-2
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
