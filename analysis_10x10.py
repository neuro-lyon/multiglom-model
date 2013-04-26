# -*- coding:utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import tables

from h5manager import get_all_attrs
from utils import to1d


plt.ion()
DB = tables.openFile('db.h5')


"""
MPS & STS

"""
IDX_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('results', '_v_attrs', 'MPS'),
             ('results', '_v_attrs', 'STS'))
IDX = get_all_attrs(DB, IDX_ATTRS)
IDX.sort()
X_IDX = list(set([i[0] for i in IDX]))
X_IDX.sort()
Y_IDX = list(set([i[1] for i in IDX]))
Y_IDX.sort()

Z_IDX = []
for i in xrange(6):
    Z_IDX.append(np.zeros((len(X_IDX), len(Y_IDX))))
for ind_rate in xrange(len(X_IDX)):
    for ind_strength in xrange(len(Y_IDX)):
        tmp_mps = IDX[to1d(ind_rate, ind_strength, len(X_IDX))][2]
        Z_IDX[0][ind_rate][ind_strength] = tmp_mps[0]
        Z_IDX[1][ind_rate][ind_strength] = tmp_mps[1]
        Z_IDX[2][ind_rate][ind_strength] = tmp_mps['whole']

        tmp_sts = IDX[to1d(ind_rate, ind_strength, len(X_IDX))][3]
        Z_IDX[3][ind_rate][ind_strength] = tmp_sts[0]
        Z_IDX[4][ind_rate][ind_strength] = tmp_sts[1]
        Z_IDX[5][ind_rate][ind_strength] = tmp_sts['whole']

# MPS plotting
IDX_FIG, IDX_AXS = plt.subplots(3)
for i in xrange(3):
    cs = IDX_AXS[i].contourf(X_IDX, Y_IDX, Z_IDX[i])
    IDX_FIG.colorbar(cs, ax=IDX_AXS[i])
plt.show()

# STS plotting
IDX_FIG, IDX_AXS = plt.subplots(3)
for i in xrange(3):
    cs = IDX_AXS[i].contourf(X_IDX, Y_IDX, Z_IDX[i + 3])
    IDX_FIG.colorbar(cs, ax=IDX_AXS[i])
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
plt.ioff()
plt.show()


DB.close()
