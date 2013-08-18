# -*- coding:utf-8 -*-

"""
Script to see if the model is like the model of [1]_ by plotting the network
frequency against the oscillation rate.

References
----------
.. [1] Fourcaud-Trocmé, N., Courtiol, E., Buonviso, N., & Voegtlin, T. (2011).
   Stability of fast oscillations in the mammalian olfactory bulb: experiments
   and modeling. Journal of physiology, Paris, 105(1-3), 59–70.
   doi:10.1016/j.jphysparis.2011.07.009

"""

import tables
import numpy as np
import matplotlib.pyplot as plt

from sys import argv
from h5manager import get_all_attrs


def plot_fI_curve(filename_data):
    data = np.loadtxt(filename_data)
    currents = data[0]
    freqs = data[1:]
    mean = freqs.mean(axis=0)
    err = freqs.std(axis=0)
    plt.errorbar(currents, mean, yerr=err)


def plot_netw_freq(db_filename, point_color, label):
    """Plot g_Ein0 against FFTMAX for the given DB."""
    db = tables.openFile(db_filename)  # Open the HDF5 database-like
    
    # Get the interesting values
    attrs_list = (('paramset', '_v_attrs', 'Input', 'g_Ein0'),
                  ('results', '_v_attrs', 'FFTMAX', 0))
    attrs = np.array(get_all_attrs(db, attrs_list))
    
    # Put them on the figure
    plt.plot(attrs[:, 0], attrs[:, 1], ' .', color=point_color, label=label)
    plt.legend(loc="upper left")

    # Finally, close the db
    db.close()


def plot_freqs(db_filename, point_color, label):
    """Plot mitral firing rate against network frequency."""
    db = tables.openFile(db_filename)

    # Get the values and arrays
    attrs_list = (('results', 'spikes_it'),
                  ('results', '_v_attrs', 'FFTMAX', 0),
                  ('paramset', '_v_attrs', 'Common'))
    attrs = get_all_attrs(db, attrs_list)
    ps_common = attrs[0][2]
    n_mitral = ps_common['N_mitral']
    simu_length = ps_common['simu_length']
    burnin = ps_common['burnin']

    # Compute the spiking rate for each simulation
    sim_values = np.ndarray((len(attrs), 2))
    for ind_simu, simu in enumerate(attrs):
        spike_times = simu[0].read()[1]
        sim_values[ind_simu][0] = get_spiking_rate(spike_times, n_mitral,
                                                   simu_length, burnin)
        sim_values[ind_simu][1] = simu[1]  # FFTMAX already computed

    # Plot the values
    plt.plot(sim_values[:, 0], sim_values[:, 1], ' .', color=point_color,
             label=label)
    plt.legend()

    # Close the DB
    db.close()


def get_spiking_rate(spike_times, n_mitral, simu_length, burnin):
    """Return the spiking rate for the whole population."""
    time_mask = (spike_times > burnin)
    return 1.*time_mask.sum()/(n_mitral*(simu_length - burnin))


def main(data_prefix):
    # Get the data
    filename_beta = data_prefix + "db40_beta_1pop_fig_netw_freq_multiproc.h5"
    filename_gamma = data_prefix + "db40_gamma_1pop_fig_netw_freq.h5"

    # fI curve
    plt.figure()
    plt.subplot(2, 1, 1)
    plot_fI_curve("data/fI_curve_data.txt")

    # Build network frequency figure
    netwfreq_axes = plt.subplot(2, 2, 1)
    plot_netw_freq(filename_beta, 'blue', u"beta (nouv. modèle)")
    plot_netw_freq(filename_gamma, 'red', u"gamma (nouv. modèle)")

    plt.xlabel("Input excitatory conductance $g_{Ein0}$ ($S m^{-2}$)")
    plt.ylabel("Network frequency $f$ (Hz)")

    # Build freq vs. freq figure
    plt.subplot(2, 2, 2, sharey=netwfreq_axes)
    plot_freqs(filename_beta, 'blue', u"beta (nouv. modèle)")
    plot_freqs(filename_gamma, 'red', u"gamma (nouv. modèle)")

    plt.xlabel("Mitral firing rate $\\nu_0$")
    plt.show()


if __name__ == '__main__':
    res = main(argv[1])
