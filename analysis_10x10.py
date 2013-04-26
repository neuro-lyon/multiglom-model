# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import tables

from h5manager import get_all_attrs
from utils import to1d

"""
TODO
- Faire tous les MPS
- Faire tous les STS

"""

plt.ion()
DB = tables.openFile('db.h5')


"""
MPS

"""
MPS_FIG, MPS_AXS = plt.subplots(3, 1)

MPS_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('results', '_v_attrs', 'MPS'))
MPS = get_all_attrs(DB, MPS_ATTRS)
MPS.sort()
X_MPS = list(set([i[0] for i in MPS]))
X_MPS.sort()
Y_MPS = list(set([i[1] for i in MPS]))
Y_MPS.sort()

Z_MPS = []
for i in xrange(3):
    Z_MPS.append(np.zeros((len(X_MPS), len(Y_MPS))))
for ind_rate in xrange(len(X_MPS)):
    for ind_strength in xrange(len(Y_MPS)):
        tmp_mps = MPS[to1d(ind_rate, ind_strength, len(X_MPS))][2]
        Z_MPS[0][ind_rate][ind_strength] = tmp_mps[0]
        Z_MPS[1][ind_rate][ind_strength] = tmp_mps[1]
        Z_MPS[2][ind_rate][ind_strength] = tmp_mps['whole']

for i in xrange(3):
    cs = MPS_AXS[i].contourf(X_MPS, Y_MPS, Z_MPS[i], MPS_AXS[i],
                            levels=np.linspace(0, 1, 11))
    MPS_FIG.colorbar(cs, ax=MPS_AXS[i])
plt.show()


"""
DEPHASAGE

"""

PHI_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('results', '_v_attrs', 'phase_angles', 0, 1))
PHI = get_all_attrs(DB, PHI_ATTRS)
PHI.sort()
X_PHI = list(set([i[0] for i in PHI]))
X_PHI.sort()
Y_PHI = list(set([i[1] for i in PHI]))
Y_PHI.sort()

Z_PHI = np.zeros((len(X_PHI), len(Y_PHI)))
for ind_rate in xrange(len(X_PHI)):
    for ind_strength in xrange(len(Y_PHI)):
        Z_PHI[ind_rate][ind_strength] = PHI[to1d(ind_rate, ind_strength, len(X_PHI))][2]

PHI_FIG = plt.figure()
PHI_CS = plt.contourf(X_PHI, Y_PHI, Z_PHI)
PHI_FIG.colorbar(PHI_CS)
plt.show()

DB.close()
