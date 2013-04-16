# -*- coding:utf-8 -*-

import tables
import numpy as np


N_SIMU = 'n_simu'


class SysState(tables.isDescription):
    run_date = tables.Time32Col()
    git_rev = tables.StringCol(40)
    comment = tables.StringCol()

class Parameters(tables.isDescription):
    Common_simu_dt = tables.FloatCol()
    Common_simu_length = tables.FloatCol()
    Common_N_subpop = tables.IntCol()
    Common_N_mitral = tables.IntCol()
    Common_inter_conn_rate = 
    Common_inter_conn_strength

    Glomerule_tau = tables.FloatCol()
    Glomerule_f = tables.FloatCol()
    Glomerule_A = tables.FloatCol()
    Glomerule_B = tables.FloatCol()
    Glomerule_C = tables.FloatCol()

    Mitral_C_m = tables.FloatCol()
    Mitral_g_L = tables.FloatCol()
    Mitral_E_L = tables.FloatCol()
    Mitral_V_r = tables.FloatCol()
    Mitral_V_t = tables.FloatCol()
    Mitral_t_refract = tables.FloatCol()

    Granule_C_m = tables.FloatCol()
    Granule_g_L = tables.FloatCol()
    Granule_E_L = tables.FloatCol()
    Granule_g_SD = tables.FloatCol()
    Granule_g_DS = tables.FloatCol()

    Input_tau_Ein = tables.FloatCol()
    Input_g_Ein0 = tables.FloatCol()
    Input_sigma_Ein = tables.FloatCol()

    Synapse_V_E = tables.FloatCol()
    Synapse_V_act_E = tables.FloatCol()
    Synapse_g_E = tables.FloatCol()
    Synapse_sigma_E = tables.FloatCol()
    Synapse_alpha_E = tables.FloatCol()
    Synapse_beta_E = tables.FloatCol()

    Synapse_V_I = tables.FloatCol()
    Synapse_V_act_I = tables.FloatCol()
    Synapse_g_I = tables.FloatCol()
    Synapse_sigma_I = tables.FloatCol()
    Synapse_alpha_I = tables.FloatCol()
    Synapse_beta_I = tables.FloatCol()


def init_data_h5(filename):
    """Initialize a data HDF5 file"""
    with tables.openFile(filename, 'w') as f:
        f.createArray('/', N_SIMU, np.array(0))


def new_simu(filename, data):
    """Put the simulation data into the HDF5 file"""
    with tables.openFile(filename, 'a') as f:
        n_simu = f.getNode('/', N_SIMU).read()
        # parse data and put them in a new group
        simu_group = f.createGroup('/', 'simu' + str(n_simu))
        # TODO change value of n_simu
