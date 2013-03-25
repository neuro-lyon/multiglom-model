from brian.stdunits import *
from brian.units import *

parameters = {
'Common':
    {'simu_dt'    : 0.05*msecond,
    'simu_length' : 2000*msecond,
    'N_subpop'    : 1,
    'N_mitral'    : 100
    },
'Glomerule':
    {'tau' : 3*msecond,
    'f'    : 2*Hz,
    'A'    : .1*second**-.5*siemens*meter**-2,
    'B'    : 10*second**-1*siemens*meter**-2,
    'C'    : 1
    },
'Mitral':
    {'C_m'      : 0.08*farad*meter**-2,
    'g_L'       : 0.87*siemens*meter**-2,
    'E_L'       : -64.5*mvolt,
    'V_r'       : -74*mvolt,
    'V_t'       : -62*mvolt,
    't_refract' : 0.2*msecond
    },
'Granule':
    {'C_m' : 0.01*farad*meter**-2,
    'g_L'  : 0.83*siemens*meter**-2,
    'E_L'  : -70*mvolt,
    'g_SD' : 1*siemens*meter**-2,
    'g_DS' : 300*siemens*meter**-2
    },
'Input':
    {'tau_Ein'  : 3*msecond,
    'g_Ein0'    : 0.6*siemens*meter**-2,
    'sigma_Ein' : .05*siemens*meter**-2*second**(-1./2)
    },
'Synapse':
    {'V_E'     : 0*mvolt,
    'V_act_E'  : 0*mvolt,
    'g_E'      : 3.5*siemens*meter**-2,
    'sigma_E'  : 0.01*mvolt,
    'alpha_E'  : 10*msecond**-1,
    'beta_E'   : 1./3*msecond**-1,

    'V_I'      : -80*mvolt,
    'V_act_I'  : -66.4*mvolt,
    'g_I'      : 10*siemens*meter**-2,
    'sigma_I'  : 0.4*mvolt,
    'alpha_I'  : 5*msecond**-1,
    'beta_I'   : 1./10*msecond**-1
    },
}
