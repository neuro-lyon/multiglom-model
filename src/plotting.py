# -*- coding:utf-8 -*-


from matplotlib import pyplot as plt, cm as cmap
from numpy import where
from brian.stdunits import *
from brian.units import *
from matplotlib.mlab import psd
from pylab import detrend_mean


def raster_plot(spikes_i, spikes_t, connection_matrix):
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
    rasterp = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
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
    rasterp.set_ylabel("Neuron number")

    # Raster histogram
    rasterhisto = plt.subplot2grid((4, 1), (3, 0), sharex=rasterp)
    nbins = spikes_t[-1] // 5e-3  # make bins of 5 ms
    rasterhisto.hist(spikes_t, bins=nbins)
    rasterhisto.set_xlabel("Time (s)")
    rasterhisto.set_ylabel("Number of spikes")

    plt.suptitle("Raster plot")

    # Connection matrix plot
    plt.figure()
    plt.imshow(connection_matrix, interpolation="nearest", extent=(0, 1, 0, 1),
               vmin=0, vmax=1)
    plt.colorbar()


def get_colorlist(n_colors, cmap_name="Paired"):
    """Get a list of n_colors colors from a matplotlib colormap."""
    colors = []
    colormap = cmap.get_cmap(cmap_name)
    assert colormap != None, cmap_name + " is not a valid colormap name."
    for i in xrange(n_colors):
        colors += [colormap(1.*i/n_colors)]
    return colors


def memb_plot_figure(monit_mt, monit_gr, rec_neurons, n_granule):
    """Membrane potentials of mitral and granule cells."""
    plt.figure()
    sub_v_mt = plt.subplot(2, 1, 1)
    for neur in rec_neurons:
        sub_v_mt.plot(monit_mt['V'].times/msecond,
                      monit_mt['V'][neur]/mvolt)
    sub_v_mt.set_xlabel('Time (ms)')
    sub_v_mt.set_ylabel('Membrane potential of mitral : V (mvolt)')

    sub_vd_gr = plt.subplot(2, 1, 2, sharex=sub_v_mt)
    for gran in xrange(n_granule):
        sub_vd_gr.plot(monit_gr['V_D'].times/msecond,
                       monit_gr['V_D'][gran]/mvolt, label="granule #" + str(gran))
    sub_vd_gr.legend()
    sub_vd_gr.set_xlabel('Time (ms)')
    sub_vd_gr.set_ylabel('Membrane potential of granule : V (mvolt)')


def granule_figure(monit_gr, pscommon):
    """Wraper to the granule figure."""
    granule_pop_figure(monit_gr['s'].values, monit_gr['s_syn_self'].values, monit_gr['s'].times, pscommon['resample_dt'], pscommon['burnin'])


def granule_pop_figure(gr_s, gr_s_syn_self, times, dt, burnin):
    """Plot a figure describing the granule activity, useful to see population synchrony."""
    plt.figure()
    n_granule = len(gr_s)

    # Granule s
    sub_s = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=3)
    for num_granule in xrange(n_granule):
        sub_s.plot(times/msecond, gr_s[num_granule],
            label="s granule #" + str(num_granule))
    sub_s.legend()
    sub_s.set_xlabel('times (ms)')
    sub_s.set_ylabel('s granule')

    # Granule s_syn_self
    sub_s_syn_self = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=3,
                                      sharex=sub_s)
    for num_granule in xrange(n_granule):
        sub_s_syn_self.plot(times/msecond, gr_s_syn_self[num_granule],
            label="s_syn_self granule #" + str(num_granule))
    sub_s_syn_self.legend()
    sub_s_syn_self.set_xlabel('times (ms)')
    sub_s_syn_self.set_ylabel('s_syn_self granule')

    # FFT max granules
    sub_fft = plt.subplot2grid((4, 4), (0, 3), rowspan=4, colspan=1)
    fft_max_freq = 200
    sig_start = where(times > burnin)[0][0]
    for num_granule in xrange(n_granule):
        power, freqs = psd(gr_s[num_granule][sig_start:], Fs=int(1/dt),
                       NFFT=int(0.5/dt), noverlap=int(0.25/dt),
                       detrend=detrend_mean)
        ind_max_freq = where(freqs <= fft_max_freq)[0][-1]
        sub_fft.plot(freqs[:ind_max_freq], power[:ind_max_freq],
                     label="FFT on granule #" + str(num_granule) + " s")
    sub_fft.legend()
    sub_fft.set_xlabel("granule s frequency (Hz)")
    sub_fft.set_ylabel('Power')


def plot_single_simulation(spikes_i, spikes_t, connection_matrix,
                           s_granule, s_syn_self, times, dt, burnin):
    """Plot figures for a single simulation"""
    # Raster plot
    raster_plot(spikes_i, spikes_t, connection_matrix)
    # Granule figure
    granule_pop_figure(s_granule, s_syn_self, times, dt, burnin)
    plt.show()
