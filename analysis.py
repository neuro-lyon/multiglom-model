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


def sts(netw_act, spikes, start, stop, keep_ratio=1./2):
    """
    Returns the STS index [1] for the given network activity.

    Parameters
    ----------
    netw_act : brian.StateMonitor.values
        Signal that represents the network activity in one variable.
    spikes : brian.SpikeMonitor
        Set of spikes during the simulation.
    start : int
        neuron index left border for the slice of neuron we want
    stop : int
        neuron index right border

    References
    ----------
    [1] Brunel & Wang, 2003

    """
    sig_size = len(netw_act)
    cut_sig = netw_act[sig_size*(1 - keep_ratio):]
    # Then compute the autocorrelation at zero time.
    autocorr = autocorr_zero(cut_sig)
    # Finally, normalize it by nu*tau
    nu = get_nspikes(spikes, keep_ratio, start, stop)/(spikes.clock.t*keep_ratio)
    return float(autocorr/(nu*TAU))  # float() to not return a Quantity object


def autocorr_zero(signal):
    """Returns the autocorrelation of the signal at zero time."""
    mean_sig = np.mean(signal)
    return np.sqrt(np.mean((signal - mean_sig)*(signal - mean_sig)))


def get_nspikes(spikes, keep_ratio, start, stop):
    """Returns the number of spikes, keeping only the last portion of the
    simulation."""
    nspikes = 0
    # Convert time to second because times in `spikes` are in second
    time_treshold = (spikes.clock.t/second)*keep_ratio
    for neur in xrange(start, stop):
        for spike_time in spikes[neur]:
            if spike_time > time_treshold:
                nspikes += 1
    return nspikes


def mps(memb_pot, start, stop, keep_ratio=1./2):
    """
    Returns the MPS index [1] of the given network.

    Parameters
    ----------
    memb_pot : StateMonitor
        Membrane potential for a whole category (eg. mitral) of neurons.
    start, stop : int
        indices of the first and last neuron to take
    keep_ratio : float
        portion of the signal to keep (right part is kept, left is dismissed)

    References
    ----------
    [1] Brunel & Wang, 2003

    """
    res = 0.
    cut_index = len(memb_pot.values[0])*(1 - keep_ratio)
    all_corr = np.corrcoef(memb_pot.values[start:stop, cut_index:])
    nneur = stop - start
    ncomb = comb(nneur, 2, exact=True)
    assert ncomb > 0, \
        "No mitral combination are possible, are you using 1 mitral?"

    for i in xrange(nneur):
        for j in xrange(i + 1, nneur):
                res += all_corr[i][j]

    return res/ncomb


def fftmax(signal, n_subpop, simu_dt, fft_max_freq=200, keep_ratio=1./2):
    """Return the peak in the FFT frequency of the signal values."""
    res = {}
    ntimes = int(len(signal.times)*(1 - keep_ratio))
    # Cut the signal values to keep_ratio
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
