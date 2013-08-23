# -*- coding:utf-8 -*-

"""
Select specific simulations from a HDF5 file and plot the granule figures and
scatter plot to compare the simulations.

"""

from plotting import get_colorlist
import matplotlib.pyplot as plt
import tables
from numpy import allclose, linspace, where
from scipy.signal import resample

import h5manager as h5m
from analysis import fftmax
from arg_parsers import ANACOMP_PARSER

def raster_plot(spikes_i, spikes_t, connection_matrix, sss):
    """Raster plot with focus on interconnection neurons.

    Parameters
    ----------
    spikes_i: array
        spike times
    spikes_t: array
        neuron number associated with spike time
    connection_matrix: array
        connection matrix of size (M mitrales, G granules)
    """
    # Raster plot
    plt.figure()
    rasterp = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    bin_connection_matrix = (connection_matrix > 0)
    n_mitral, n_subpop = connection_matrix.shape
    n_mitral_per_subpop = n_mitral/n_subpop

    # Make a mapping, neuron: {spike times}
    spike_map = {}
    for neur, time in zip(spikes_i, spikes_t):
        if spike_map.has_key(neur):
            spike_map[int(neur)].append(time)
        else:
            spike_map[int(neur)] = [time]

    # Plotting
    colors = get_colorlist(n_subpop)
    for ind_subpop in xrange(n_subpop):
        subpop_start = ind_subpop*n_mitral_per_subpop
        subpop_stop  = subpop_start + n_mitral_per_subpop
        subpop_color = colors[ind_subpop]
        downline = subpop_start
        upline   = subpop_stop - 1
        for ind_neuron in xrange(subpop_start, subpop_stop):
            neur_connections = bin_connection_matrix[ind_neuron]
            # Getting the neuron spike times, if it spiked
            if spike_map.has_key(ind_neuron):
                spikes = spike_map[ind_neuron]
            else:
                spikes = []
            # Plotting the spikes for that neuron
            if neur_connections.sum() > 1:  # if the neuron is connected to more than one granule
                dark_color = [i/1.5 for i in subpop_color[:-1]]
                rasterp.plot(spikes, [upline]*len(spikes), ' .',
                         color=dark_color, mew=0)
                upline -= 1
            else:
                rasterp.plot(spikes, [downline]*len(spikes), ' .',
                         color=subpop_color, mew=0)
                downline += 1

    # Some plotting enhancement
    margin = 0.01
    if len(spikes_t) > 0:
        spikes_t_last = spikes_t[-1]
    else:
        spikes_t_last = 0.
    x_overplot = margin*spikes_t_last
    y_overplot = margin*n_mitral
    rasterp.set_xlim((-x_overplot), (spikes_t_last + x_overplot))
    rasterp.set_ylim(-y_overplot, n_mitral + y_overplot)
    rasterp.set_ylabel("Indice des neurones")

    # Raster histogram
    rasterhisto = plt.subplot2grid((5, 1), (3, 0), sharex=rasterp)
    nbins = spikes_t[-1] // 5e-3  # make bins of 5 ms
    rasterhisto.hist(spikes_t, bins=nbins)
    rasterhisto.set_ylabel("Nombre de spikes")

    # S syn self
    raster_sss = plt.subplot2grid((5, 1), (4, 0), sharex=rasterp)
    plt.plot(linspace(spikes_t[0], spikes_t[-1], sss.T.shape[0]), sss.T)
    raster_sss.set_xlabel("Temps (s)")
    raster_sss.set_ylabel("s_syn_self")


if __name__ == '__main__':

    ARGS = ANACOMP_PARSER.parse_args()

    DB_FILENAME = ARGS.data_file
    DB = tables.openFile(DB_FILENAME)

    PLOT_MEMB_POT = ARGS.plot_mp

    # Get all simulation data
    ATTRS = [('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'simu_length'),
             ('paramset', '_v_attrs', 'Common', 'simu_dt'),
             ('results', 's_granule'),
             ('results', 's_syn_self'),
             ('results', '_v_attrs'),
             ('results', 'spikes_it'),
             ('paramset', 'arrays', 'mtgr_connections'),
    ]
    if PLOT_MEMB_POT:
        ATTRS.append(('results', 'mean_memb_pot'))

    ALL_SIMU_ATTRS = h5m.get_all_attrs(DB, ATTRS)
    ALL_SIMU_ATTRS.sort()

    SET_INTERCO_RATE = list(set(s[0] for s in ALL_SIMU_ATTRS))
    SET_INTERCO_RATE.sort()
    SET_INTERCO_STRENGTH = list(set(s[1] for s in ALL_SIMU_ATTRS))
    SET_INTERCO_STRENGTH.sort()


    # Some useful functions
    def deq(a, b, delta):
        """Delta equality"""
        return abs(a - b) < delta

    class SignalRepack:
        def __init__(self, values, times):
            self.values = values
            self.times = times

    # Get the specific simulations
    def get_interco(simu, interco_rate, interco_strength):
        return deq(simu[0], interco_rate, 0.02) and deq(simu[1], interco_strength, 0.02)

    BURNIN = 1.
    SELECTED_RATES = [22]
    SELECTED_STRENGTH = [6]
    REDO_FFTMAX = []
    for rate in SELECTED_RATES:
        for strength in SELECTED_STRENGTH:
            interco_rate = SET_INTERCO_RATE[rate]
            interco_strength = SET_INTERCO_STRENGTH[strength]
            filtered_simu = filter(lambda x: allclose(x[0], interco_rate) \
                                             and allclose(x[1], interco_strength),
                                   ALL_SIMU_ATTRS)
            if len(filtered_simu) > 1:
                print "Warning: more than one simulation satisfy your criteria,",
                print "selecting the first one."
            simu = filtered_simu[0]

            # Granule plot
            gr_s = simu[4].read()
            gr_s_syn_self = simu[5].read()
            simu_length = float(simu[2])
            resample_dt = simu_length/len(gr_s[1])
            times = linspace(0., simu_length, len(gr_s[0]))
            sig_start = where(times > BURNIN)[0][0]
            mtgr_connections = simu[8].read()
            # granule_pop_figure(gr_s, gr_s_syn_self, times, resample_dt, BURNIN)

            # Raster plot
            spikes_it = simu[7].read()
            raster_plot(spikes_it[0], spikes_it[1], mtgr_connections, gr_s_syn_self)

            # # Membrane potential
            # if PLOT_MEMB_POT:
            #     memb_potentials = simu[9].read()
            #     labels = ("Mean. non-interco.", "Mean. interco")
            #     plt.figure()
            #     for ind_mp, memb_pot in enumerate(memb_potentials):
            #         label = labels[ind_mp % 2]
            #         label += " glom #" + str(ind_mp/2)
            #         plt.plot(times, memb_pot, label=label)
            #         plt.legend()
            #     plt.xlabel("Time (s)")
            #     plt.ylabel("Membrane potential (V)")

            # FFT max peak
            signal = SignalRepack(gr_s, times)
            REDO_FFTMAX.append(fftmax(signal, 2, resample_dt, sig_start))

            print 'rate:', rate, 'strength:', strength, REDO_FFTMAX
    plt.show()

    # Finally close the DB
    DB.close()
