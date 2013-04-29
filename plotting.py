# -*- coding:utf-8 -*-

from matplotlib import pyplot as plt, cm as cmap
from scipy.fftpack import fft, fftfreq
from brian.stdunits import *
from brian.units import *


def raster_plot(spike_monitor, n_subpop):
    """Raster plot with a different color for each sub-population."""
    n_neurons = len(spike_monitor.spiketimes)
    n_neuron_per_subpop = n_neurons/n_subpop
    colors = get_colorlist(n_subpop)
    plt.figure()
    for subpop in xrange(n_subpop):
        subpop_color = colors[subpop]
        for neur in xrange(n_neuron_per_subpop):
            abs_neuron = neur + subpop*n_neuron_per_subpop
            spikes = spike_monitor[abs_neuron]
            plt.plot(spikes/msecond, [abs_neuron]*len(spikes), ' .',
                     color=subpop_color, mew=0)
    margin = 0.01
    x_overplot = margin*spike_monitor.clock.end
    y_overplot = margin*n_neurons
    plt.xlim((-x_overplot)/msecond,
             (spike_monitor.clock.end + x_overplot)/msecond)
    plt.ylim(-y_overplot, n_neurons + y_overplot)
    plt.xlabel("time (ms)")
    plt.ylabel("neuron number")


def get_colorlist(n_colors, cmap_name="gist_ncar"):
    """Get a list of `n_colors` color from a matplotlib colormap."""
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
    """Parameters to/from the granule, useful to see population synchrony."""
    for gr in xrange(pscommon['N_subpop']):
        granule_pop_figure(monit_gr['s'].values, monit_gr['s_syn_self'].values, gr, monit_gr['s'].times, pscommon['simu_dt'])


def granule_pop_figure(gr_s, gr_s_syn_self, num_granule, times, dt):
    """Plot a figure describing the granule activity."""
    plt.figure()

    # Granule s & s_syn activities
    sub_s = plt.subplot(1, 2, 1)
    sub_s.plot(times/msecond, gr_s_syn_self[num_granule],
        label="s_syn_self granule #" + str(num_granule))
    sub_s.plot(times/msecond, gr_s[num_granule],
        label="s granule #" + str(num_granule))
    sub_s.legend()
    sub_s.set_xlabel('time (ms)')
    sub_s.set_ylabel('s & s_syn_self granule #' + str(num_granule))

    # Granule FFT
    sub_synchro = plt.subplot(1, 2, 2)
    fft_max_freq = 200
    ntimes = len(times)
    freqs = fftfreq(ntimes, dt)
    fft_max_freq_index = next(f for f in xrange(len(freqs)) if freqs[f] > fft_max_freq)
    fft_sig = abs(fft(gr_s[num_granule] - (gr_s[num_granule]).mean())[:fft_max_freq_index])

    sub_synchro.plot(freqs[:fft_max_freq_index], fft_sig)
    sub_synchro.set_xlabel("granule #"+str(num_granule)+" 's' frequency (Hz)")
    sub_synchro.set_ylabel('Power')
