.. _workflow:

Workflow
========

How to run a single simulation
------------------------------

1. Create or edit file that has all the parameters for the model.
   See template files in ``src/paramsets``.
2. Run ``python multiglom_network.py path_to/parameter_file.py``.

How to run multiple simulations
-------------------------------

1. Create all parameter files.
   To do so you need to edit the bottom of ``src/utils.py``.
2. Commit your changes with git.
3. Run the simulations with :mod:`run_multip.py`
4. Use :meth:`collect_h5_to_db` to get all the simulations together.
