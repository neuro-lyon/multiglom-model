# -*- coding:utf-8 -*-

"""
Analysis
========

Provides different functions to analyse the neuronal network activity.

Synchronization
---------------
- Spike Train Synchrony (STS) index
- Membrane Potential Synchrony (MPS) index
- FFT peak

"""
import numpy as np
from model import PARAMETERS as ps

from scipy.misc import comb
from scipy.fftpack import fft, fftfreq
from scipy.signal import resample
from scipy import argmax

from brian.stdunits import *
from brian.units import *

if ps:
    PSIN = ps['Input']
    TAU  = PSIN['tau_Ein']


def sts(netw_act, spikes, neur_start, neur_stop, sig_start, time_start):
    """
    Returns the STS index [1] for the given network activity.

    Parameters
    ----------
    netw_act: brian.StateMonitor.values
        signal that represents the network activity in one variable
    spikes: brian.SpikeMonitor
        set of spikes during the simulation
    neur_start: int
        neuron index left border for the slice of neuron we want
    neur_stop: int
        neuron index right border

    References
    ----------
    [1] Brunel & Wang, 2003

    """
    cut_sig = netw_act[sig_start:]
    # Then compute the autocorrelation at zero time.
    autocorr = autocorr_zero(cut_sig)
    # Finally, normalize it by nu*tau
    nspikes = get_nspikes(spikes, time_start, neur_start, neur_stop)
    nu = nspikes/(spikes.clock.end - time_start)
    return float(autocorr/(nu*TAU))  # float() to not return a Quantity object


def autocorr_zero(signal):
    """Returns the autocorrelation of the signal at zero time."""
    mean_sig = np.mean(signal)
    return np.sqrt(np.mean((signal - mean_sig)*(signal - mean_sig)))


def get_nspikes(spikes, time_treshold, start, stop):
    """Returns the number of spikes, keeping only the last portion of the
    simulation."""
    nspikes = 0
    for neur in xrange(start, stop):
        for spike_time in spikes[neur]:
            if spike_time > time_treshold:
                nspikes += 1
    return nspikes


def mps(memb_pot, start, stop, sig_start):
    """
    Returns the MPS index [1] of the given network.

    Parameters
    ----------
    memb_pot : StateMonitor
        Membrane potential for a whole category (eg. mitral) of neurons.
    start, stop : int
        indices of the first and last neuron to take

    References
    ----------
    [1] Brunel & Wang, 2003

    """
    res = 0.
    all_corr = np.corrcoef(memb_pot.values[start:stop, sig_start:])
    nneur = stop - start
    ncomb = comb(nneur, 2, exact=True)
    assert ncomb > 0, \
        "No mitral combination are possible, are you using 1 mitral?"

    for i in xrange(nneur):
        for j in xrange(i + 1, nneur):
                res += all_corr[i][j]

    return res/ncomb


def fftmax(signal, n_subpop, simu_dt, sig_start, fft_max_freq=200):
    """Return the peak in the FFT frequency of the signal values."""
    res = {}
    ntimes = int(len(signal.times[sig_start:]))
    cut_signal = signal.values[:, ntimes:]

    freqs = fftfreq(ntimes, simu_dt)
    fft_max_freq_index = next(f for f in xrange(len(freqs)) if freqs[f] > fft_max_freq)

    # Compute FFT for each subpopulation
    for unit in xrange(n_subpop):
        fft_sig = abs(fft(cut_signal[unit]-(cut_signal[unit]).mean())[:fft_max_freq_index])
        ind_max_freq = argmax(fft_sig)
        res[unit] = freqs[ind_max_freq]

    # Compute FFT for the whole population by the mean of activities
    mean_signal = np.mean(cut_signal, axis=0)
    fft_sig = abs(fft(mean_signal - mean_signal.mean())[:fft_max_freq_index])
    ind_max_freq = argmax(fft_sig)
    res['mean'] = freqs[ind_max_freq]

    return res


def crosscorr_phase_angle(sig1, sig2, x, max_length=10000):
    """Return the cross correlation phase angle between 2 signals

    Parameters
    ----------
    sig1 : array
        signal of length L
    sig2 : array
        another signal of length L
    x : array
        time axis for the signals sig1 and sig2
    max_length : int
        maximum length for the signals, signals are resampled otherwise
    """
    assert len(sig1) == len(sig2) == len(x), \
        "The signals don't have the same length."
    sig_length = len(sig1)
    # Resample if signal is too big thus slowing down correlation computation
    if sig_length > max_length:
        sig1, x = resample(sig1, max_length, x)
        sig2 = resample(sig2, max_length)
        sig_length = max_length
    corr = np.correlate(sig1, sig2, mode="same")
    xmean = sig_length/2
    return float(argmax(corr) - xmean)/sig_length*x[-1]  # *x[-1] to scale


def peak_dist_index(sig1, sig2, xaxis=None):
    """Return the mean and std of peak-distances between the signals"""
    peak_dist = get_dist(sig1, sig2, xaxis)
    return np.mean(peak_dist), np.std(peak_dist)


def peak_dist_circ_index(sig1, sig2):
    """Return the *circular* mean and std of the distances between peaks"""
    # Make the distances directional
    first_sig, _ = get_ordered_sig([sig1, sig2])
    n_peaks = len(get_ind_local_max(first_sig))
    mean_peak_dist = len(first_sig)/n_peaks
    peak_dist = np.array(get_dist(sig1, sig2))
    peak_dist = peak_dist*2.*np.pi/mean_peak_dist
    # Apply circular statistcs to the directional distances
    return circ_mean(peak_dist), circ_disp(peak_dist)


def circ_disp(sig):
    """Statistical circular dispertion

    References
    ----------
    [1] http://cran.r-project.org/web/packages/CircStats/
    """
    n = len(sig)
    scos = np.sum(np.cos(sig))
    ssin = np.sum(np.sin(sig))
    root = np.sqrt(scos*scos + ssin*ssin)
    rbar = root/n
    var  = 1 - rbar
    return var


def circ_mean(sig):
    """Statistical circular mean

    References
    ----------
    [1] http://cran.r-project.org/web/packages/CircStats/
    """
    sinr = np.sum(np.sin(sig))
    cosr = np.sum(np.cos(sig))
    return np.arctan2(sinr, cosr)


def get_dist(sig1, sig2, xaxis=None):
    """Return the distances between the peaks of two signals"""
    distances = []
    max_sig1 = get_ind_local_max(sig1)
    max_sig2 = get_ind_local_max(sig2)
    first_sig, second_sig = get_ordered_sig((max_sig1, max_sig2))
    ind_peak_fs = 0
    ind_peak_ss = 0
    last_peak = False

    while ind_peak_fs < len(first_sig) and ind_peak_ss < len(second_sig):
        # Take the indexes of the peaks
        peak_fs = first_sig[ind_peak_fs]
        peak_ss = second_sig[ind_peak_ss]
        if ind_peak_fs == len(first_sig) - 1:  # if this is the last peak,
            peak_fs_next = sig1[-1]  # make up a last peak at the very end
            last_peak = True
        else:
            peak_fs_next = first_sig[ind_peak_fs + 1]
        # If an x-axis is given, get the peak x values from it
        if xaxis is not None:
            peak_fs = xaxis[peak_fs]
            peak_fs_next = xaxis[peak_fs_next]
            peak_ss = xaxis[peak_ss]

        dist_intra_fs = peak_fs_next - peak_fs
        dist_inter = peak_ss - peak_fs
        # No SS peak in between two FS peaks
        if dist_intra_fs < dist_inter and not last_peak:
            ind_peak_fs += 1
        # There is one or more SS peaks in between two FS peaks, or it's the last peak
        else:
            dist_left = peak_ss - peak_fs
            dist_right = peak_ss - peak_fs_next
            if abs(dist_left) < abs(dist_right):
                distances.append(dist_left)
            else:
                distances.append(dist_right)
            ind_peak_ss += 1

    return distances


def get_ind_local_max(sig, treshold_ratio=0.25):
    """Return indexes of the local maxima of sig.

    Parameters
    ----------
    sig : list or np.array
        signal to get the local maxima from
    treshold_ratio : float
        ratio relative to the distance median to keep the peaks,
        e.g. if the distance median is 10 and the ratio is 0.25 then only
        peak distances above 2.5 will be kept.
    """
    peak_indexes =  np.nonzero((np.diff(sig[:-1]) > 0) & (np.diff(sig[1:]) < 0))[0] + 1
    peak_distances = np.diff(peak_indexes)
    peak_dist_median = np.median(peak_distances)
    suptreshold_peaks = (peak_distances > treshold_ratio*peak_dist_median)
    return peak_indexes[suptreshold_peaks]


def get_ordered_sig(sig_list):
    """Return the signals ordered according to their first value"""
    first_values = []
    
    # Get the first values of each signal to sort them
    for i, sig in enumerate(sig_list):
        first_values.append((sig[0], i))
    first_values.sort()

    ordered_list = []
    for fv in first_values:
        sig = sig_list[fv[1]]
        ordered_list.append(sig)

    return ordered_list


def slope(x1, y1, x2, y2):
    """Return the slope between (x1, y1) and (x2, y2)."""
    return (y2 - y1)/(x2 - x1)


def sign(val):
    """Return the sign of the value: -1 if negative, +1 if positive or null."""
    if val < 0:
        return -1
    elif val >= 0:
        return 1
