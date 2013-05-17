# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import tables

from h5manager import get_all_attrs
from utils import to1d


FILENAME='data/db30x30_two_glom_beta_new_ps_interco_strength0_1_interco_rate0_1.h5'
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
    fig.subplots_adjust(bottom=0.25)
    for ind_subplot, ind_data in enumerate(data_range):
        cs = axs[ind_subplot].imshow(data[ind_data],
                                     origin="lower",
                                     norm=data_norm,
                                     interpolation="nearest",
                                     extent=axes_extent,
                                     aspect="auto")
        axs[ind_subplot].set_title(fig_title)
    cb_axs = fig.add_axes([0.125, 0.1, 0.9 - 0.125, 0.03])
    fig.colorbar(cs, cax=cb_axs, orientation="horizontal")


"""
Common Attributes

"""
COMMON_ATTRS = (('paramset', '_v_attrs', 'Common'),)
COMMON = get_all_attrs(DB, COMMON_ATTRS)
COMMON = COMMON[0][0]
SIMU_LENGTH = COMMON['simu_length']
SIMU_DT = COMMON['simu_dt']
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
            ('results', 'spikes_it'))
SR = get_all_attrs(DB, SR_ATTRS)
SR.sort()
X_SR = list(set(i[0] for i in SR))
X_SR.sort()
Y_SR = list(set(i[1] for i in SR))
Y_SR.sort()

START_TIME = SIMU_LENGTH/2.
Z_SR = []
SPIKES_IT = np.array([])
for i in xrange(3):
    Z_SR.append(np.zeros((len(X_SR), len(Y_SR))))
for ind_rate in xrange(len(X_SR)):
    for ind_strength in xrange(len(Y_SR)):
        SPIKES_IT = np.array(SR[to1d(ind_rate, ind_strength, len(X_SR))][2])
        nspikes = []
        for pop in xrange(N_SUBPOP):
            neuron_range=[pop*N_MITRAL_PER_SUBPOP, (pop + 1)*N_MITRAL_PER_SUBPOP]
            valid_spikes=(SPIKES_IT[0,:]>=neuron_range[0])&(SPIKES_IT[0,:]<neuron_range[1])
            valid_times = (SPIKES_IT[1,valid_spikes] > float(START_TIME))
            nspikes.append(valid_times.sum()/(N_MITRAL_PER_SUBPOP*(SIMU_LENGTH - START_TIME)))
        nspikes_whole = 1.*(SPIKES_IT[1,:] >= START_TIME).sum()/(N_MITRAL * (SIMU_LENGTH - START_TIME))
        nspikes.append(nspikes_whole)
        Z_SR[0][ind_rate][ind_strength] = nspikes[0]
        Z_SR[1][ind_rate][ind_strength] = nspikes[1]
        Z_SR[2][ind_rate][ind_strength] = nspikes[2]

# Plotting
SR_FIG, SR_AXS = plt.subplots(1, 3, figsize=(9, 3))
SR_NORM = colors.normalize(np.amin(Z_SR), np.amax(Z_SR))
plot_run(SR_FIG, "Spiking Rate", SR_AXS, Z_SR, range(3), SR_NORM, IMSHOW_EXTENT)


"""
Closing DB and finally plotting

"""
DB.close()
plt.show()
