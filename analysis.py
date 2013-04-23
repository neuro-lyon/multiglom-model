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
from scipy import argmax

from brian.stdunits import *
from brian.units import *

PSIN = ps['Input']
TAU  = PSIN['tau_Ein']


def sts(netw_act, spikes, start, stop, keep_ratio=3./4):
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


def mps(memb_pot, start, stop):
    """
    Returns the MPS index [1] of the given network.

    Parameters
    ----------
    memb_pot : StateMonitor
        Membrane potential for a whole category (eg. mitral) of neurons.

    References
    ----------
    [1] Brunel & Wang, 2003

    """
    res = 0.
    all_corr = np.corrcoef(memb_pot.values[start:stop])
    nneur = stop - start
    ncomb = comb(nneur, 2, exact=True)
    assert ncomb > 0, \
        "No mitral combination are possible, are you using 1 mitral?"
    for i in xrange(nneur):
        for j in xrange(nneur):
            if j > i:
                res += all_corr[i][j]
    return res/ncomb


def fftmax(signal, n_subpop, simu_dt, fft_max_freq=200):
    """Return the peak in the FFT frequency of the signal values."""
    res = []
    ntimes = len(signal.times)
    freqs = fftfreq(ntimes, simu_dt)
    fft_max_freq_index = next(f for f in xrange(len(freqs)) if freqs[f] > fft_max_freq)

    for unit in xrange(n_subpop):
        fft_sig = abs(fft(signal[unit]-(signal[unit]).mean())[:fft_max_freq_index])
        ind_max_freq = argmax(fft_sig)
        res.append(freqs[ind_max_freq])

    return res
