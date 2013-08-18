#!/usr/bin/env python2.7
# -*- coding:utf-8 -*-
"""
This is a script to get the fI curve for the mitral cells.

"""
import brian_no_units
from brian import *


def main(args):
    import model_utils as mutils
    # Set the parameters from the specified file BEFORE any model.* import
    import model
    mutils.set_model_ps(args.psfile)

    import numpy as np
    from scipy.signal import resample

    from model.mitral_cells import MitralCells

    # Reset old stuff from Brian memory
    clear(erase=True, all=True)
    defaultclock.reinit()

    # Initialize random generator (necessary mainly for parallel simulations)
    np.random.seed()
    
    """
    Parameters
    ----------
    """
    psmt     = model.PARAMETERS['Mitral']
    pscommon = model.PARAMETERS['Common']

    n_trials = 10
    n_mitral = pscommon['N_mitral']*n_trials

    defaultclock.dt = pscommon['simu_dt']
    simu_length     = pscommon['simu_length']


    """
    Population Initialization
    -------------------------
    """

    # Mitral cells
    mt = MitralCells()
    mt_supp_eqs =  {'var': ['+I_electrode'],
                    'eqs': ['I_electrode : amp*meter**-2']}
    mt.add_eqs(supp_eqs=mt_supp_eqs)

    mt.make_pop(n_mitral)

    # Inject some currents for the fI curve
    i_start = 0.*amp*meter**-2
    i_stop = 0.01*amp*meter**-2
    currents = np.linspace(i_start, i_stop, n_mitral//n_trials)
    for ind_trial in xrange(n_trials):
        pop_start = n_mitral//n_trials * ind_trial
        pop_stop = pop_start + n_mitral//n_trials
        mt.pop.I_electrode[pop_start:pop_stop] = currents

    mt.pop.V = (psmt['V_t'] - psmt['V_r'])*np.random.random_sample(np.shape(mt.pop.V)) \
               + psmt['V_r']


    """
    Simulation Monitoring
    ---------------------
    """

    monit_spikes = SpikeMonitor(mt.pop, record=True)

    """
    Running Simulation
    ------------------
    """
    # Gathering simulation objects
    netw = Network(mt.pop, monit_spikes)

    # Simulation run
    if args.no_brian_output:
        report_output = None
    else:
        report_output = "text"
    netw.run(simu_length, report=report_output)


    """
    Simulation output
    -----------------
    """
    spiketimes = monit_spikes.spiketimes
    freqs = []
    for neuron in spiketimes:
        freqs += [len(spiketimes[neuron])/pscommon['simu_length']]
    # Reshaping, one line by trial
    freqs_by_trial = np.reshape(freqs, (n_trials, n_mitral//n_trials))

    res = np.ndarray((n_trials + 1, n_mitral//n_trials))
    res[0] = currents[:n_mitral//n_trials]
    res[1:] = freqs_by_trial
    np.savetxt("fI_curve_data.txt", res)
    return res

if __name__ == '__main__':
    # Argument parsing
    from arg_parsers import SIM_PARSER
    args = SIM_PARSER.parse_args()

    # Run script
    res = main(args)
