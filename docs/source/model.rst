Model
=====

The multi-glomeruli model is defined by 4 main classes : :class:`Glomerule`,
:class:`MitralCells`, :class:`GranuleCells` and :class:`Synapse`.
Each one of these is mainly a set of equations describing a number of neurons
(except for :class:`Synapse` which is just a set of equations).

To make one of :class:`Glomerule`, :class:`MitralCells` or
:class:`GranuleCells`, you must :

1. initialize it by making an instance of the class
2. setting the equation model by calling :meth:`set_eqs`
3. making a :class:`NeuronGroup` by calling :meth:`make_pop`

Note that at this point, the populations are not connected to each other.

Glomerule
---------
.. currentmodule:: glomerule
.. autoclass:: Glomerule
   :members:

   .. automethod:: __init__

Mitral Cells
------------
.. currentmodule:: mitral_cells
.. autoclass:: MitralCells
   :members:

   .. automethod:: __init__

Granule Cells
-------------
.. currentmodule:: granule_cells
.. autoclass:: GranuleCells
   :members:

   .. automethod:: __init__

Synapse
-------
.. currentmodule:: synapse
.. autoclass:: Synapse
   :members:

   .. automethod:: __init__
