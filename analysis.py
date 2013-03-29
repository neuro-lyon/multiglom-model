#!/usr/bin/env python2
# -*- coding:utf-8 -*-

"""
Analysis
========

Provides different functions to analyse the neuronal network activity.

Synchronization
---------------
- Spike Train Synchrony (STS) index
- Membrane Potential Synchrony (MPS) index
- Cross-correlation
- Gabor Transform

"""
import numpy as np
import model.parameters as ps

from scipy.signal import resample
from scipy.misc import comb

from brian.stdunits import *
from brian.units import *

PSIN = ps.Input()
TAU  = PSIN.tau_Ein
N_SAMPLES = 1000


def sts(netw_act, spikes, nsamples=N_SAMPLES):
    """
    Returns the STS index [1] for the given network activity.

    Parameters
    ----------
    netw_act : brian.StateMonitor
        Signal that represents the network activity in one variable.
    spikes : brian.SpikeMonitor
        Set of spikes during the simulation.
    nsamples : int
        Number of samples when resampling the signal.

    References
    ----------
    [1] Brunel & Wang, 2003

    """
    keep_ratio = 3./4
    sig_size = len(netw_act[0])
    cut_sig = netw_act[0][sig_size*(1 - keep_ratio):]
    # Fist, resample the signal
    signal = resample(cut_sig, nsamples)
    # Then compute the autocorrelation at zero time.
    autocorr = autocorr_zero(signal)
    # Finally, normalize it by sqrt(nu)*tau
    nu = get_nspikes(spikes, keep_ratio)/(spikes.clock.t*keep_ratio)
    return autocorr/(nu*TAU)


def autocorr_zero(signal):
    """Returns the autocorrelation of the signal at zero time."""
    mean_sig = np.mean(signal)
    return np.sqrt(np.mean((signal - mean_sig)*(signal - mean_sig)))


def get_nspikes(spikes, keep_ratio):
    """Returns the number of spikes, keeping only the last portion of the
    simulation."""
    spike_times = spikes.it[1]
    start_time = 0
    while spike_times[start_time]*second < (1 - keep_ratio)*spikes.clock.t:
        start_time += 1
    return int(len(spike_times[start_time:]))



def mps(memb_pot):
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
    all_corr = np.corrcoef(memb_pot.values)
    nneur = len(memb_pot.record)
    for i in xrange(nneur):
        for j in xrange(nneur):
            if j > i:
                res += all_corr[i][j]
    return res/comb(nneur, 2, exact=True)
