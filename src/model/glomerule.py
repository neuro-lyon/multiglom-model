# -*- coding:utf-8 -*-

from brian import *

from model import PARAMETERS as ps

# If the model parameters have been initialized
if ps != None:
    PSIN = ps['Input']
    tau_Ein   = PSIN['tau_Ein']
    g_Ein0    = PSIN['g_Ein0']
    sigma_Ein = PSIN['sigma_Ein']

    PSINOSC = ps['InputOscillation']
    f    = PSINOSC['f']
    C    = PSINOSC['C']


class Glomerule:
    """Ensemble of glomeruli, driving input to mitral cells.

    A :class:`Glomerule` reprensents one glomerule, but it is made of multiple
    smaller units. Each one of these smaller units is connected down to a mitral
    cell.

    The real role of a :class:`Glomerule` is acting as network input.

    Attributes
    ----------
    eqs_model : Brian.Equations
        Equations that define the granule model
    pop : Brian.NeuronGroup
        Object that represents a set of neurons
    """

    def __init__(self):
        """Create an empty glomerule.

        Notes
        -----
        :meth:`add_eqs` and :meth:`make_pop` must be called after initializing
        the granule.
        """
        self.eqs_model = Equations()
        self.pop = None
    
    def add_eqs(self, oscillating=False):
        """Add standard equations to the glomerule model.

        Parameters
        ----------
        oscillating : bool, optional
            True to generate oscillating input
        """
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
        """Make the NeuronGroup for a population of glomeruli.

        Parameters
        ----------
        N : int
            Total number of small glomerule units in the whole simulation.
            N = nb. glomeruli * nb. mitral cell per glomerular column

        Notes
        -----
        This must be called after all the other construction methods.
        """
        self.pop = NeuronGroup(N, model=self.eqs_model)
