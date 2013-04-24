from brian.stdunits import *
from brian.units import *
PARAMETERS = {'Mitral': {'V_t': -62.0 * mvolt, 'V_r': -74.0 * mvolt, 't_refract': 0.2 * msecond, 'C_m': 0.08 * metre ** -4 * kilogram ** -1 * second ** 4 * amp ** 2, 'g_L': 0.87 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2, 'E_L': -64.5 * mvolt}, 'Synapse': {'V_E': 0.0 * volt, 'sigma_I': 0.4 * mvolt, 'beta_E': 0.333333333333 * khertz, 'beta_I': 0.1 * khertz, 'V_I': -80.0 * mvolt, 'sigma_E': 10.0 * uvolt, 'V_act_I': -66.4 * mvolt, 'alpha_E': 10.0 * khertz, 'V_act_E': 0.0 * volt, 'alpha_I': 5.0 * khertz, 'g_I': 10.0 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2, 'g_E': 3.5 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2}, 'InputOscillation': {'C': 1, 'f': 2.0 * hertz}, 'Granule': {'g_L': 0.83 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2, 'g_SD': 1.0 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2, 'C_m': 0.01 * metre ** -4 * kilogram ** -1 * second ** 4 * amp ** 2, 'g_DS': 300.0 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2, 'E_L': -70.0 * mvolt}, 'Common': {'simu_length': 3.0 * second, 'N_subpop': 2, 'inter_conn_strength': {0: {1: 0.0}, 1: {0: 0.0}}, 'N_mitral': 100, 'simu_dt': 50.0 * usecond, 'inter_conn_rate': {0: {1: 0.77777777777777768}, 1: {0: 0.77777777777777768}}}, 'Input': {'tau_Ein': 3.0 * msecond, 'sigma_Ein': 0.05 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2, 'g_Ein0': 1.0 * metre ** -4 * kilogram ** -1 * second ** 3 * amp ** 2}}
