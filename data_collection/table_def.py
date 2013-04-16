# -*- coding:utf-8 -*-
"""
Describe the tables in the main HDF5 file.

"""


from tables import Time32Col, StringCol, IntCol, FloatCol, isDescription


class SysState(isDescription):
    """System state"""
    run_date = Time32Col()
    git_rev = StringCol(40)
    comment = StringCol()
    numpy_seed = IntCol()


class Parameters(isDescription):
    """Parameter set"""
    Common_simu_dt = FloatCol()
    Common_simu_length = FloatCol()
    Common_N_subpop = IntCol()
    Common_N_mitral = IntCol()

    Glomerule_tau = FloatCol()
    Glomerule_f = FloatCol()
    Glomerule_A = FloatCol()
    Glomerule_B = FloatCol()
    Glomerule_C = FloatCol()

    Mitral_C_m = FloatCol()
    Mitral_g_L = FloatCol()
    Mitral_E_L = FloatCol()
    Mitral_V_r = FloatCol()
    Mitral_V_t = FloatCol()
    Mitral_t_refract = FloatCol()

    Granule_C_m = FloatCol()
    Granule_g_L = FloatCol()
    Granule_E_L = FloatCol()
    Granule_g_SD = FloatCol()
    Granule_g_DS = FloatCol()

    Input_tau_Ein = FloatCol()
    Input_g_Ein0 = FloatCol()
    Input_sigma_Ein = FloatCol()

    Synapse_V_E = FloatCol()
    Synapse_V_act_E = FloatCol()
    Synapse_g_E = FloatCol()
    Synapse_sigma_E = FloatCol()
    Synapse_alpha_E = FloatCol()
    Synapse_beta_E = FloatCol()

    Synapse_V_I = FloatCol()
    Synapse_V_act_I = FloatCol()
    Synapse_g_I = FloatCol()
    Synapse_sigma_I = FloatCol()
    Synapse_alpha_I = FloatCol()
    Synapse_beta_I = FloatCol()
