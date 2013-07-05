# -*- coding:utf-8 -*-
"""
Model related utilities.

"""

import importlib
import model
import numpy as np

from brian import StateMonitor, SpikeMonitor
from utils import path_to_modline


def set_model_ps(filepath, dicname='PARAMETERS'):
    """Set the model parameters given the parameter set.

    Parameters
    ----------
    filepath : str
        File with the parameter set
    dicname : str, optional
        Variable in :attr:`filepath` that is the dictionnary of parameters
    """
    psmod = importlib.import_module(path_to_modline(filepath))
    model.PARAMETERS = getattr(psmod, dicname)


def intrapop_connections(n_mitral, n_granule, n_subpop, n_mitral_per_subpop):
    """Build a connection matrix for intra sub-population connections.

    Parameters
    ----------
    n_mitral : int
        Number of mitral cells in the network
    n_granule : int
        Number of granule cells in the network
    n_subpop : int
        Number of sub-population in the network
    n_mitral_per_subpop : int
        Number of mitral cells per sub-population in the network

    Returns
    -------
    np.ndarray
        Connection matrix between mitral cells (rows) and granule cells (cols)
    """
    resmat = np.zeros((n_mitral, n_granule))
    for i_subpop in xrange(n_subpop):
        start = i_subpop*n_mitral_per_subpop
        stop  = start + n_mitral_per_subpop
        resmat[start:stop, i_subpop] = 1.
    return resmat


def interpop_connections(mat_connections, n_mitral, n_subpop, n_mitral_per_subpop, inter_conn_rate, inter_conn_strength, homeostasy=False):
    """
    Adds inter sub-population connections.

    Parameters
    ----------
    mat_connections : np.ndarray
        Connection matrix between mitral cells and granule cells
    n_mitral : int
        Number of mitral cells in the network
    n_subpop : int
        Number of sub-population in the network
    n_mitral_per_subpop : int
        Number of mitral cells per sub-population in the network
    inter_conn_rate : float
        Interconnection rate [0, 1]
    inter_conn_strength : float
        Interconnection strength [0, 1]
    homeostasy : bool, optional
        True to set homeostatic connections.
        Default is False.

    Returns
    -------
    np.ndarray
        Connection matrix between mitral cells (rows) and granule cells (cols)
    """
    if homeostasy:
        init_total = 1.*mat_connections.sum(axis=0)
        tr_init_total = 1.*mat_connections.sum(axis=1)

    res_mat = mat_connections
    n_granule = n_subpop
    subpop_start = np.zeros((n_subpop)) # which is the first non interconnected mitral from each subpop
    for mtpop in inter_conn_rate:
        assert mtpop >= 0 and mtpop < n_subpop, \
            "Incorrect mitral sub-population number "+str(mtpop)+" for inter-connectivity."
        for grpop in inter_conn_rate[mtpop]:
            if grpop != mtpop:
                assert grpop >= 0 and grpop < n_granule, \
                    "Incorrect granule sub-population number "+str(grpop)+" for inter-connectivity."
                conn = inter_conn_rate[mtpop][grpop]
                assert conn >= 0 and conn <= 1, "Connectivity must be in [0, 1]."
                nlinks = int(n_mitral_per_subpop*conn)
                newconn = np.zeros((n_mitral_per_subpop, 1))
                for i in xrange(nlinks):
                    try:
                        newconn[i + subpop_start[mtpop]] = inter_conn_strength[mtpop][grpop]
                    except Exception, e:
                        print "Likely, too much connections to have no overlap, rewrite the code !"
                        print "EXCEPTION:", e
                        exit(1)
                subpop_start[mtpop] += nlinks
                start = mtpop*n_mitral_per_subpop
                stop  = start + n_mitral_per_subpop
                res_mat[start:stop, grpop] = newconn[:, 0]
    
    if not(homeostasy):
        return res_mat, res_mat.transpose()
    else:
        # Connection strengthes are normalized such that each neuron (mitral or granule) receive the same total amount of excitation or inhibition 
        # (i.e. the same total of synaptic conductance) as if they were no interglomerular connections
        tr_res_mat = res_mat.transpose().copy()
        tr_res_mat_norm = tr_init_total/tr_res_mat.sum(axis=0)  # here the numerator is 1. because there is only one granule
        tr_res_mat = tr_res_mat*tr_res_mat_norm
        res_mat_norm = init_total/res_mat.sum(axis=0)
        res_mat = res_mat*res_mat_norm
        return res_mat, tr_res_mat


def monit(pop, params, timestep, reclist=True, spikes=False):
    """Returns a dictionnary of monitors for a population.

    Parameters
    ----------
    pop : brian.NeuronGroup
        Neuron population to record
    params : list
        List of parameters to record
    timestep : int
        Timestep of the simulatoin
    reclist : bool or list, optional
        List of neurons to record, or True to record them all.
        Default is True
    spikes : bool
        True to put a brian.SpikeMonitor to record, instead of a StateMonitor

    Returns
    -------
    dict
        Dictionnary of parameter name (keys) and monitors (values)
    """
    res = {}
    for pname in params:
        res[pname] = StateMonitor(pop, pname, record=reclist, timestep=timestep)
    if spikes:
        res['spikes'] = SpikeMonitor(pop, record=True)
    return res
