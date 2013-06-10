#!/usr/bin/env python2.7
# -*- coding:utf-8 -*-
"""
Multi-glomerular Network Simulation
===================================

This is a script that runs a simulation of a small OB model.
It consists of an input signal and some OB columns. A column is a set of one
glomerule, its connected mitral cells, and the granule cells.

Script Overview
---------------
1. Import the model and necessary stuff (like numpy)
2. Get the set of parameters (see inside the paramsets directory)
3. Initialize the different cell populations:
    - glomeruli
    - (synapses between granule and mitral)
    - mitral cells
    - granule cells
4. Connects the different cell populations:
    - glomeruli and mitral cells
    - mitral cells and granule cells
5. Set some monitors on the simulation
6. Run the simulation
7. Output simulation information and indexes.
8. Plots

"""
import brian_no_units
from brian import *


def main(args):
    import model_utils as mutils
    # Set the parameters from the specified file BEFORE any model.* import
    import model
    mutils.set_model_ps(args.psfile)

    import numpy as np
    import analysis
    import plotting
    from utils import print_dict, pairs
    from scipy.signal import resample

    from model.glomerule import Glomerule
    from model.mitral_cells import MitralCells
    from model.synapse import Synapse
    from model.granule_cells import GranuleCells

    # Reset old stuff from Brian memory
    clear(erase=True, all=True)
    defaultclock.reinit()

    # Initialize random generator (necessary mainly for parallel simulations)
    np.random.seed()
    
    """
    Parameters
    ----------
    Get the parameter values from the `ps` module, which in turn gets the values
    from the file specified in parameters.py.

    Set some aliases for the different cell population sizes.
    Also check that there is an even number of cells for each column.

    Finally set some simulation parameters.

    """
    psmt     = model.PARAMETERS['Mitral']
    psgr     = model.PARAMETERS['Granule']
    pscommon = model.PARAMETERS['Common']

    n_mitral    = pscommon['N_mitral']
    n_glomeruli = n_granule = n_subpop = pscommon['N_subpop']

    # check to have an even number of mitral in each sub-population
    assert n_mitral % n_subpop == 0, \
           "N_mitral is not a multiple of the number of sub-populations N_subpop."
    n_mitral_per_subpop = n_mitral/n_subpop

    defaultclock.dt = pscommon['simu_dt']
    simu_length     = pscommon['simu_length']


    """
    Population Initialization
    -------------------------
    1. glomeruli
    *. synapses between granule and mitral cells
    3. mitral cells
    4. granule cells

    """
    # Glomeruli
    glom = Glomerule()
    glom.add_eqs()
    glom.make_pop(n_glomeruli*n_mitral_per_subpop)

    # Synapses (granule -- mitral)
    synexc = Synapse(synapse_type='exc') # excitatory synapse
    synexc.set_eqs_model()

    syninhib = Synapse(synapse_type='inhib') # inhibitory synapse
    syninhib.set_eqs_model()

    # Mitral cells
    mt = MitralCells()
    mt_supp_eqs =  {'var': ['- I_syn', '- g_input*V'],
                    'eqs': [synexc.get_eqs_model(),
                            Equations("g_input : siemens*meter**-2")]}
    mt.add_eqs(supp_eqs=mt_supp_eqs)
    mt.make_pop(n_mitral)
    mt.pop.V = (psmt['V_t'] - psmt['V_r'])*np.random.random_sample(np.shape(mt.pop.V)) \
               + psmt['V_r']

    # Granule Cells
    gr = GranuleCells()
    gr_supp_eqs = {'var': ['-I_syn'],
                   'eqs': [syninhib.get_eqs_model()]}
    gr.add_eqs(supp_eqs=gr_supp_eqs)
    gr.make_pop(n_granule)
    gr.pop.V_D = psgr['E_L']
    gr.pop.V_S = psgr['E_L']


    """
    Connecting Populations
    ----------------------
    1. Glomeruli and mitral cells 
    2. Mitral cells and granule cells

    """
    # Connecting mitral cells to glomeruli
    glmt_connections = diag(ones(n_mitral))

    # Glomeruli--Mitral interactions
    @network_operation(when='start')
    def mt_input():
        mt.pop.g_input = dot(glom.pop.g, glmt_connections)

    # Connecting sub-population of mitral cells to granule cells
    mtgr_connections = mutils.intrapop_connections(n_mitral, n_granule, n_subpop, n_mitral_per_subpop)

    # Inter subpopulation connectivities
    inter_conn_rate = pscommon['inter_conn_rate']
    inter_conn_strength = pscommon['inter_conn_strength']
    homeostasy = pscommon['homeostasy']
    mtgr_connections, grmt_connections = mutils.interpop_connections(mtgr_connections, n_mitral, n_subpop,
                                            n_mitral_per_subpop, inter_conn_rate, inter_conn_strength,homeostasy)
    # Mitral--Granule interactions
    @network_operation(when='start')
    def graded_synapse():
        """Computes granule and mitral s_syn"""
        mt.pop.state('T')[:] = 0.
        mt.pop.state('T')[mt.pop.get_refractory_indices()] = 1.
        gr.pop.s_syn = dot(mt.pop.s, mtgr_connections)
        mt.pop.s_syn = dot(gr.pop.s, grmt_connections)

    @network_operation(when='start')
    def sum_s():
        """Computes granule self s_syn (for its glomerular column only)"""
        for subpop in xrange(n_subpop):
            start = subpop*n_mitral_per_subpop
            stop  = start + n_mitral_per_subpop
            gr.pop.s_syn_self[subpop] = sum(mt.pop.state('s')[start:stop])

    @network_operation(when='after_groups')
    def keep_reset():
        mt.pop.state('V')[mt.pop.get_refractory_indices()] = psmt['V_r']


    """
    Simulation Monitoring
    ---------------------
    Monitor state variables for the different populations.

    """
    glom_ps = ('g')
    mt_ps   = ('s', 's_syn', 'V')
    gr_ps   = ('V_D', 's_syn', 's', 's_syn_self')

    # Simulation monitors
    rec_neurons = True  # Must be set to True if we want accurate MPS and STS
    monit_glom = mutils.monit(glom.pop, glom_ps, reclist=rec_neurons)
    monit_mt   = mutils.monit(mt.pop, mt_ps, reclist=rec_neurons, spikes=True)
    monit_gr   = mutils.monit(gr.pop, gr_ps)


    """
    Running Simulation
    ------------------
    Create Network object and put everything simulation related in it.
    Then run this network.

    """
    # Gathering simulation objects
    netw = Network(glom.pop, mt.pop, gr.pop,
                   mt_input, graded_synapse, keep_reset, sum_s,
                   [m for m in monit_glom.values()],
                   [m for m in monit_mt.values()],
                   [m for m in monit_gr.values()])

    # Simulation run
    if args.no_brian_output:
        report_output = None
    else:
        report_output = "text"
    netw.run(simu_length, report=report_output)


    """
    Information Output
    ------------------

    """
    if args.full_ps:
        print 'Full set of parameters:'
        print_dict(model.PARAMETERS)

    burnin = pscommon['burnin']
    times = monit_gr['s'].times
    sig_start = where(times > burnin)[0][0]

    sts_indexes = {}
    mps_indexes = {}
    fftmax = {}
    mps_indexes['whole'] = analysis.mps(monit_mt['V'], 0, n_mitral, sig_start)
    gr_s_syn_whole = np.zeros(monit_gr['s_syn'][0].shape)


    # MPS and STS computation for subpopulation
    for subpop in xrange(n_subpop):
        start = subpop*n_mitral_per_subpop
        stop  = start + n_mitral_per_subpop
        sts = analysis.sts(monit_gr['s_syn'][subpop], monit_mt['spikes'], start, stop, sig_start, burnin)
        sts_indexes[subpop] = sts
        gr_s_syn_whole += monit_gr['s_syn'][subpop]
        mps = analysis.mps(monit_mt['V'], start, stop, sig_start)
        mps_indexes[subpop] = mps

    # STS for the whole population
    sts_whole_activity = np.zeros(monit_gr['s_syn'][0].shape)
    for subpop in xrange(n_subpop):
        sts_whole_activity += monit_gr['s_syn'][subpop]  # TODO pas tres sur s'il faut ajouter ou concaténer
    sts_indexes['whole'] = analysis.sts(sts_whole_activity, monit_mt['spikes'], 0, n_mitral, sig_start, burnin)

    # FFT Max index
    fftmax = analysis.fftmax(monit_gr['s'], n_subpop, pscommon['simu_dt'], sig_start)

    # Peak distances index
    peak_distances = {}
    if n_subpop > 1:
        for sub_i, sub_j in pairs(n_subpop):
            sig1 = monit_gr['s_syn_self'][sub_i]
            sig2 = monit_gr['s_syn_self'][sub_j]
            if not peak_distances.has_key(sub_i):
                peak_distances[sub_i] = {}
            pd_index = analysis.peak_dist_circ_index(sig1, sig2)
            peak_distances[sub_i][sub_j] = {}
            peak_distances[sub_i][sub_j]['mean'] = pd_index[0]
            peak_distances[sub_i][sub_j]['disp'] = pd_index[1]

    if not args.no_summary:
        print '\nParameters: using', args.psfile

        print 'Populations:', n_subpop, 'glomerular columns;',
        print n_mitral, 'mitral cells;', n_granule, 'granule cells.'

        print 'Times:', simu_length, 'of simulation; dt =', defaultclock.dt, '.'

        print 'Indexes: STS =', sts_indexes, '\nMPS =', mps_indexes
        print 'FFT peaks (Hz):', fftmax
        print 'Peak distances index:', peak_distances

    """
    Plotting
    --------
    Plot monitored variables and a scatter plot.

    """
    if not args.no_plot:
        # Raster plot
        spikes_it = monit_mt['spikes'].it
        plotting.raster_plot(spikes_it[0], spikes_it[1], mtgr_connections)
        # Membrane potentials
        if not rec_neurons:  # if we only have a couple of recorded neurons
            plotting.memb_plot_figure(monit_mt, monit_gr, rec_neurons, n_granule)
        # Granule synapses
        plotting.granule_figure(monit_gr, pscommon)
        show()


    """
    Simulation records
    ------------------

    Put numpy arrays in var `results` to save them into the simulation record.
    Note: the variable must be monitored by Brian.

    """
    # Add parameters
    ps_arrays = {'mtgr_connections': (mtgr_connections,
                            "Connection matrix from mitral (rows) to granules (columns)")}

    # Add results
    array_spikes_it = np.array((monit_mt['spikes'].it[0],
                                monit_mt['spikes'].it[1]))
    results = {}

    # Mean inputs
    mean_inputs = np.ndarray((n_glomeruli, monit_glom['g'].values.shape[1]))
    for glom in xrange(n_glomeruli):
        start_subpop = glom*n_mitral_per_subpop
        stop_subpop = start_subpop + n_mitral_per_subpop
        mean_inputs[glom] = np.mean(monit_glom['g'].values[start_subpop:stop_subpop], axis=0)

    # Mean membrane potentials
    mean_memb_pot = np.ndarray((n_glomeruli*2, monit_mt['V'].values.shape[1]))
    bin_interco_matrix = (mtgr_connections > 0.)
    interco_neurons = (bin_interco_matrix.sum(axis=1) > 1)
    for glom in xrange(n_glomeruli):
        start_subpop = glom*n_mitral_per_subpop
        stop_subpop = start_subpop + n_mitral_per_subpop
        # Get subpopulation membrane potentials and interconnected neurons
        subpop_memb_pot = monit_mt['V'].values[start_subpop:stop_subpop]
        subpop_interco_neurons = interco_neurons[start_subpop:stop_subpop]
        # Compute one mean for interconnected neurons and another for the other neurons
        mean_pop = np.mean(subpop_memb_pot[~subpop_interco_neurons], axis=0)
        mean_pop_interco = np.mean(subpop_memb_pot[subpop_interco_neurons], axis=0)
        mean_memb_pot[glom*2] = mean_pop
        mean_memb_pot[glom*2 + 1] = mean_pop_interco

    results['data'] = {'spikes_it': [array_spikes_it,
                           "Spikes: one array for the neuron number, another one for the spike times."],
                       'input': [mean_inputs,
                           "Mean network input conductance value for each glomerule."],
                       's_granule': [monit_gr['s'].values,
                           "Variable 's' of the granules."],
                       's_syn_self': [monit_gr['s_syn_self'].values,
                           "Variable 's_syn' for the granule, without  integrating the mitral 's' from other subpopulations."],
                       'mean_memb_pot': [mean_memb_pot,
                            "Mean membrane potential. For each subpop: one mean for the interconnected neurons and one mean for the non-interconnected neurons."]}

    # Resample some result signals
    new_dt = 5e-4  # in second
    do_not_resample = ('spikes_it')  # list of array names to NOT resample
    for signal in results['data']:
        if signal not in do_not_resample:
            raw_signal = results['data'][signal][0]
            n_points = simu_length/new_dt
            resampled_signal = resample(raw_signal, n_points, axis=1)
            results['data'][signal][0] = resampled_signal


    results['indexes'] = {'MPS': mps_indexes, 'STS': sts_indexes, 'FFTMAX': fftmax,
                          'peak_distances': peak_distances}

    return {'set': model.PARAMETERS, 'arrays': ps_arrays}, results


if __name__ == '__main__':
    # Argument parsing
    from arg_parsers import SIM_PARSER
    args = SIM_PARSER.parse_args()

    # Run script
    ps, res = main(args)
