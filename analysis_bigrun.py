# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import tables

from h5manager import get_all_attrs
from utils import to1d
from arg_parsers import ANABIGRUN_PARSER

ARGS = ANABIGRUN_PARSER.parse_args()
FILENAME = ARGS.data_file
DB = tables.openFile(FILENAME)


"""
Utility Functions

"""
def get_all_min_max(array_list):
    """Return the min and max between all arrays."""
    big_arr = np.concatenate(array_list)
    return (np.amin(big_arr), np.amax(big_arr))


def plot_run(fig, fig_title, axs, data, data_range, data_norm, axes_extent):
    """Plot a figure with 3 subplots (subpop 1, 2, pop)"""
    fig.suptitle(fig_title)
    fig.subplots_adjust(bottom=0.25)
    for ind_subplot, ind_data in enumerate(data_range):
        cs = axs[ind_subplot].imshow(data[ind_data],
                                     origin="lower",
                                     norm=data_norm,
                                     interpolation="nearest",
                                     extent=axes_extent,
                                     aspect="auto")
    cb_axs = fig.add_axes([0.125, 0.1, 0.9 - 0.125, 0.03])
    fig.colorbar(cs, cax=cb_axs, orientation="horizontal")


"""
Common Attributes

"""
COMMON_ATTRS = (('paramset', '_v_attrs', 'Common'),)
COMMON = get_all_attrs(DB, COMMON_ATTRS)
COMMON = COMMON[0][0]
SIMU_LENGTH = float(COMMON['simu_length'])
N_MITRAL = COMMON['N_mitral']
N_SUBPOP = COMMON['N_subpop']
N_MITRAL_PER_SUBPOP = N_MITRAL/N_SUBPOP


"""
MPS & STS

"""
IDX_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('results', '_v_attrs', 'MPS'),
             ('results', '_v_attrs', 'STS'))
IDX = get_all_attrs(DB, IDX_ATTRS)
IDX.sort()
X_IDX = list(set([i[0] for i in IDX]))
X_IDX.sort()
Y_IDX = list(set([i[1] for i in IDX]))
Y_IDX.sort()

IMSHOW_EXTENT = (X_IDX[0], X_IDX[-1], Y_IDX[0], Y_IDX[-1])

Z_IDX = []
for i in xrange(6):
    Z_IDX.append(np.zeros((len(X_IDX), len(Y_IDX))))
for ind_rate in xrange(len(X_IDX)):
    for ind_strength in xrange(len(Y_IDX)):
        tmp_mps = IDX[to1d(ind_rate, ind_strength, len(X_IDX))][2]
        Z_IDX[0][ind_rate][ind_strength] = tmp_mps[0]
        Z_IDX[1][ind_rate][ind_strength] = tmp_mps[1]
        Z_IDX[2][ind_rate][ind_strength] = tmp_mps['whole']

        tmp_sts = IDX[to1d(ind_rate, ind_strength, len(X_IDX))][3]
        Z_IDX[3][ind_rate][ind_strength] = tmp_sts[0]
        Z_IDX[4][ind_rate][ind_strength] = tmp_sts[1]
        Z_IDX[5][ind_rate][ind_strength] = tmp_sts['whole']

# MPS plotting
MPS_FIG, MPS_AXS = plt.subplots(1, 3, figsize=(9, 3))
MPS_MIN_MAX = get_all_min_max([Z_IDX[i] for i in xrange(3)])
MPS_NORM = colors.normalize(MPS_MIN_MAX[0], MPS_MIN_MAX[1])
plot_run(MPS_FIG, "MPS", MPS_AXS, Z_IDX, range(3), MPS_NORM, IMSHOW_EXTENT)

# STS plotting
STS_FIG, STS_AXS = plt.subplots(1, 3, figsize=(9, 3))
STS_MIN_MAX = get_all_min_max([Z_IDX[i] for i in xrange(3, 6)])
STS_NORM = colors.normalize(STS_MIN_MAX[0], STS_MIN_MAX[1])
plot_run(STS_FIG, "STS", STS_AXS, Z_IDX, range(3, 6), STS_NORM, IMSHOW_EXTENT)


"""
FFT MAX

"""
FFT_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('results', '_v_attrs', 'FFTMAX'))
FFT = get_all_attrs(DB, FFT_ATTRS)
FFT.sort()
X_FFT = list(set(i[0] for i in FFT))
X_FFT.sort()
Y_FFT = list(set(i[1] for i in FFT))
Y_FFT.sort()

Z_FFT = []
for i in xrange(3):
    Z_FFT.append(np.zeros((len(X_FFT), len(Y_FFT))))
for ind_rate in xrange(len(X_FFT)):
    for ind_strength in xrange(len(Y_FFT)):
        tmp_fft = FFT[to1d(ind_rate, ind_strength, len(X_FFT))][2]
        Z_FFT[0][ind_rate][ind_strength] = tmp_fft[0]
        Z_FFT[1][ind_rate][ind_strength] = tmp_fft[1]
        if type(tmp_fft) == type({}):  # The first two big runs didn't record FFT mean, so skip it
            Z_FFT[2][ind_rate][ind_strength] = tmp_fft['mean']

# Plotting
FFT_FIG, FFT_AXS = plt.subplots(1, 3, figsize=(9, 3))
FFT_NORM = colors.normalize(np.amin(Z_FFT), np.amax(Z_FFT))
plot_run(FFT_FIG, "FFT", FFT_AXS, Z_FFT, range(3), FFT_NORM, IMSHOW_EXTENT)


"""
Spiking Rate

"""
SR_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
            ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
            ('results', 'spikes_it'),
            ('paramset', 'arrays', 'mtgr_connections'),)
SR = get_all_attrs(DB, SR_ATTRS)
SR.sort()
X_SR = list(set(i[0] for i in SR))
X_SR.sort()
Y_SR = list(set(i[1] for i in SR))
Y_SR.sort()

START_TIME = 1.
Z_SR = []
Z_SR_DIFF = []
SPIKES_IT = np.array([])
for i in xrange(3):
    Z_SR.append(np.zeros((len(X_SR), len(Y_SR))))
for i in xrange(6):
    Z_SR_DIFF.append(np.zeros((len(X_SR), len(Y_SR))))

for ind_rate in xrange(len(X_SR)):
    for ind_strength in xrange(len(Y_SR)):
        ind_simu = to1d(ind_rate, ind_strength, len(X_SR))
        SPIKES_IT = np.array(SR[ind_simu][2])
        mtgr_connections = SR[ind_simu][3].read()
        bin_connection_matrix = (mtgr_connections > 0.)
        nspikes = []

        for pop in xrange(N_SUBPOP):
            neuron_range = [pop*N_MITRAL_PER_SUBPOP, (pop + 1)*N_MITRAL_PER_SUBPOP]
            valid_spikes = (SPIKES_IT[0, :] >= neuron_range[0]) & (SPIKES_IT[0, :] < neuron_range[1])
            valid_times = (SPIKES_IT[1, valid_spikes] > float(START_TIME))
            nspikes.append(valid_times.sum()/(N_MITRAL_PER_SUBPOP*(SIMU_LENGTH - START_TIME)))

            # interco. vs non-interco. neurons
            subpop_interco_neurons = []
            subpop_non_interco_neurons = []
            for ind_neur in xrange(neuron_range[0], neuron_range[1]):
                # Check if the neuron is connected to more than one granule
                if bin_connection_matrix[ind_neur].sum() > 1:
                    subpop_interco_neurons.append(ind_neur)
                else:
                    subpop_non_interco_neurons.append(ind_neur)
            # Get the spikes of the neurons
            time_mask = SPIKES_IT[1] >= START_TIME
            interco_mask = np.in1d(SPIKES_IT[0], subpop_interco_neurons)
            non_interco_mask = np.in1d(SPIKES_IT[0], subpop_non_interco_neurons)
            nspikes_interco = (time_mask & interco_mask).sum()
            nspikes_non_interco = (time_mask & non_interco_mask).sum()
            # Derive the spiking rates
            n_mitrals_interco = len(subpop_interco_neurons)
            n_mitrals_non_interco = len(subpop_non_interco_neurons)
            time_window = (SIMU_LENGTH - START_TIME)
            if n_mitrals_interco != 0:
                interco_sr = nspikes_interco/(n_mitrals_interco*time_window)
            else:
                interco_sr = 0
            if n_mitrals_non_interco != 0:
                non_interco_sr = nspikes_non_interco/(n_mitrals_non_interco*time_window)
            else:
                non_interco_sr = 0
            Z_SR_DIFF[pop][ind_rate][ind_strength] = interco_sr
            Z_SR_DIFF[pop + N_SUBPOP + 1][ind_rate][ind_strength] = non_interco_sr

        nspikes_whole = 1.*(SPIKES_IT[1,:] >= START_TIME).sum()/(N_MITRAL * (SIMU_LENGTH - START_TIME))
        nspikes.append(nspikes_whole)

        Z_SR[0][ind_rate][ind_strength] = nspikes[0]
        Z_SR[1][ind_rate][ind_strength] = nspikes[1]
        Z_SR[2][ind_rate][ind_strength] = nspikes[2]

        # Compute spiking rate for all pop and split interco. and non-interco. neurons
        all_interco_mask = (bin_connection_matrix.sum(axis=1) > 1)
        all_interco_neur = np.where(all_interco_mask)[0]
        interco_mask = np.in1d(SPIKES_IT[0], all_interco_neur)

        all_non_interco_neur = np.where(~all_interco_mask)[0]
        non_interco_mask = np.in1d(SPIKES_IT[0], all_non_interco_neur)

        time_mask = SPIKES_IT[1] >= START_TIME
        nspikes_interco = (time_mask & interco_mask).sum()
        nspikes_non_interco = (time_mask & non_interco_mask).sum()

        if len(all_interco_neur) != 0:
            interco_sr = nspikes_interco/(len(all_interco_neur)*(SIMU_LENGTH - START_TIME))
        else:
            interco_sr = 0
        if len(all_non_interco_neur) != 0:
            non_interco_sr = nspikes_non_interco/(len(all_non_interco_neur)*(SIMU_LENGTH - START_TIME))
        else:
            non_interco_sr = 0

        Z_SR_DIFF[N_SUBPOP][ind_rate][ind_strength] = interco_sr
        Z_SR_DIFF[2*N_SUBPOP + 1][ind_rate][ind_strength] = non_interco_sr

# Plotting
SR_FIG, SR_AXS = plt.subplots(1, 3, figsize=(9, 3))
SR_NORM = colors.normalize(np.amin(Z_SR), np.amax(Z_SR))
plot_run(SR_FIG, "Spiking Rate", SR_AXS, Z_SR, range(3), SR_NORM, IMSHOW_EXTENT)
# Plotting Spiking Rate for interco. and non-interco. neurons
SR_DIFF_FIG, SR_DIFF_AXS = plt.subplots(1, 6)
SR_DIFF_NORM = colors.normalize(np.amin(Z_SR_DIFF), np.amax(Z_SR_DIFF))
plot_run(SR_DIFF_FIG, "Spiking Rate (interco. vs non-interco.)", SR_DIFF_AXS,
         Z_SR_DIFF, range(6), SR_DIFF_NORM, IMSHOW_EXTENT)

"""
Peak distances

"""
# Assuming only 2 subpops
PD_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
            ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
            ('results', '_v_attrs', 'peak_distances'))
PD = get_all_attrs(DB, PD_ATTRS)
PD.sort()
X_PD = list(set(i[0] for i in PD))
X_PD.sort()
Y_PD = list(set(i[1] for i in PD))
Y_PD.sort()

START_TIME = SIMU_LENGTH/2.
Z_PD = []
for i in xrange(2):
    Z_PD.append(np.zeros((len(X_PD), len(Y_PD))))
for ind_rate in xrange(len(X_PD)):
    for ind_strength in xrange(len(Y_PD)):
        tmp_pd = PD[to1d(ind_rate, ind_strength, len(X_PD))][2][0][1]
        Z_PD[0][ind_rate][ind_strength] = tmp_pd['mean']
        Z_PD[1][ind_rate][ind_strength] = tmp_pd['disp']

# Plotting
PD_FIG, PD_AXS = plt.subplots(1, 2, figsize=(6, 3))
for ind_subplot, ind_data in enumerate(range(2)):
    if ind_subplot == 0:  # if it's the mean, which is circular
        color = "hsv"  # use the hsv colormap, which is circular
        pd_norm = colors.normalize(-np.pi, np.pi)
    else:
        color = None
        pd_norm = colors.normalize(np.amin(Z_PD[1]), np.amax(Z_PD[1]))
    cs = PD_AXS[ind_subplot].imshow(Z_PD[ind_data],
                                 origin="lower",
                                 interpolation="nearest",
                                 extent=IMSHOW_EXTENT,
                                 aspect="auto",
                                 norm=pd_norm,
                                 cmap=color)
    if ind_subplot == 0:
        index_type = "mean"
    elif ind_subplot == 1:
        index_type = "disp"
    PD_AXS[ind_subplot].set_title("Peak Distances index (%s)" % index_type )
    PD_FIG.colorbar(cs, ax=PD_AXS[ind_subplot], orientation="horizontal")

"""
Closing DB and finally plotting

"""
DB.close()
plt.show()
