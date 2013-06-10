# -*- coding:utf-8 -*-

import tables
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

from h5manager import get_all_attrs

FILENAME_BETA = "data/db40_beta_1pop_fig_netw_freq.h5"
FILENAME_GAMMA = "data/db40_gamma_1pop_fig_netw_freq.h5"


def plot_netw_freq(db_filename, point_color, label):
    """Plot `g_Ein0` against `FFTMAX` for the given DB."""
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
    sim_values = np.ndarray((len(attrs), len(attrs[0])))
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
    return time_mask.sum()/(n_mitral*simu_length)


# Build network frequency figure
FIG_NETW_FREQ = plt.figure()

plot_netw_freq(FILENAME_BETA, 'blue', "beta")
plot_netw_freq(FILENAME_GAMMA, 'red', "gamma")

plt.xlabel("Input excitatory conductance $g_{Ein0}$ (S $m^{-2}$)")
plt.ylabel("Network frequency $f$ (Hz)")


# Build freq vs. freq figure
FIG_FREQS = plt.figure()

plot_freqs(FILENAME_BETA, 'blue', "beta")
plot_freqs(FILENAME_GAMMA, 'red', "gamma")

plt.xlabel("Mitral firing rate $\\nu_0$")
plt.ylabel("Network frequency $f$ (Hz)")
plt.show()
