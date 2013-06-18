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
    """Set the model parameters given the parameter set in `filepath`."""
    psmod = importlib.import_module(path_to_modline(filepath))
    model.PARAMETERS = getattr(psmod, dicname)


def intrapop_connections(n_mitral, n_granule, n_subpop, n_mitral_per_subpop):
    """Connection matrix for intra sub-population connections."""
    resmat = np.zeros((n_mitral, n_granule))
    for i_subpop in xrange(n_subpop):
        start = i_subpop*n_mitral_per_subpop
        stop  = start + n_mitral_per_subpop
        resmat[start:stop, i_subpop] = 1.
    return resmat


def interpop_connections(mat_connections, n_mitral, n_subpop, n_mitral_per_subpop, inter_conn_rate, inter_conn_strength, homeostasy=False):
    """
    Adds inter sub-population connections.

    Changes on the connection matrix are done in place.

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
            assert grpop >= 0 and grpop < n_granule, \
                "Incorrect granule sub-population number "+str(grpop)+" for inter-connectivity."
            conn = inter_conn_rate[mtpop][grpop]
            assert conn >= 0 and conn <= 1, "Connectivity must be in [0, 1]."
            nlinks = int(n_mitral_per_subpop*conn)
            newconn = np.zeros((n_mitral_per_subpop, 1))
            for i in xrange(nlinks):
                try:
                    newconn[i + subpop_start[mtpop]] = inter_conn_strength[mtpop][grpop]
                except:
                    print "Likely, too much connections to have no overlap, rewrite the code !"
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
    """Returns a dictionnary of monitors for the population."""
    res = {}
    for pname in params:
        res[pname] = StateMonitor(pop, pname, record=reclist, timestep=timestep)
    if spikes:
        res['spikes'] = SpikeMonitor(pop, record=True)
    return res
