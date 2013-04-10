# -*- coding:utf-8 -*-

from pylab import *
from brian.stdunits import *
from brian.units import *

def raster_plot(spike_monitor, n_subpop):
    """Raster plot different color for each subpop."""
    pass

def get_colorlist(n_colors, cmap_name):
    colors = []
    colormap = get_cmap(cmap_name)
    for i in xrange(n_colors):
        colors += [colormap(1.*i/n_colors)]
    return colors

def memb_plot_figure(monit_mt, monit_gr, rec_neurons, n_granule):
    figure()
    sub_v_mt = subplot(2, 1, 1)
    for neur in rec_neurons:
        sub_v_mt.plot(monit_mt['V'].times/msecond,
                      monit_mt['V'][neur]/mvolt)
    sub_v_mt.set_xlabel('Time (ms)')
    sub_v_mt.set_ylabel('Membrane potential of mitral : V (mvolt)')

    sub_vd_gr = subplot(2, 1, 2, sharex=sub_v_mt)
    for gran in xrange(n_granule):
        sub_vd_gr.plot(monit_gr['V_D'].times/msecond,
                       monit_gr['V_D'][gran]/mvolt, label="granule #" + str(gran))
    sub_vd_gr.legend()
    sub_vd_gr.set_xlabel('Time (ms)')
    sub_vd_gr.set_ylabel('Membrane potential of granule : V (mvolt)')
