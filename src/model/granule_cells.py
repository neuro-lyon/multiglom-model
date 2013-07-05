# -*- coding:utf-8 -*-

from brian import Equations, NeuronGroup
from brian.units import *
from brian.stdunits import *

from model import PARAMETERS as ps

# If the model parameters have been initialized
if ps != None:
    PSGR = ps['Granule']
    C_m  = PSGR['C_m']
    g_L  = PSGR['g_L']
    E_L  = PSGR['E_L']
    g_SD = PSGR['g_SD']
    g_DS = PSGR['g_DS']

class GranuleCells:
    """A single entity acting as a population of granule cells.

    This class represents a population of granule cells combined into a single
    :class:`Granule`. Each glomerular column has its own granule.

    The granule dynamics are defined according to a set of equations.

    Attributes
    ----------
    eqs_model : Brian.Equations
        Equations that define the granule model
    pop : Brian.NeuronGroup
        Object that represents a set of neurons
    """

    def __init__(self):
        """Create an empty population of granule cells.

        Notes
        -----
        :meth:`add_eqs` and :meth:`make_pop` must be called after initializing
        the granule.
        """
        self.eqs_model = Equations()
        self.pop = None

    def add_eqs(self, supp_eqs=None):
        """Add the standard equations of the granule model.

        Parameters
        ----------
        supp_eqs : Brian.Equations, optional
            Supplementary equations to the model

        See Also
        --------
        MitralCells.add_eqs
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
        """Makes the NeuronGroup for a population of N granule cells.

        Parameters
        ----------
        N : int
            Number of granule in the whole simulation.

        Notes
        -----
        This must be called after all the other construction methods.
        """
        self.pop = NeuronGroup(N, model=self.eqs_model, implicit=True)

