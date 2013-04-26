# -*- coding:utf-8 -*-
from pylab import *
from h5manager import *
from utils import to1d
import tables

"""
TODO
- Faire tous les MPS
- Faire tous les STS

"""

DB = tables.openFile('db.h5')

# Get data
DATA = []
for g in DB.walkGroups():
    try:
        pset = g.paramset._v_attrs
        res = g.results._v_attrs
    except tables.NoSuchNodeError:
        pass
    else:
        new_data_item = {}

        # Interconnection rate and strength
        common = pset['Common']
        interco_strength = common['inter_conn_strength'][0][1]
        interco_rate = common['inter_conn_rate'][0][1]
        new_data_item['interco'] = {}
        new_data_item['interco']['strength'] = interco_strength
        new_data_item['interco']['rate'] = interco_rate

        # MPS
        new_data_item['MPS'] = res['MPS']

        # MPS
        new_data_item['STS'] = res['STS']

        # Add data item
        DATA.append(new_data_item)



# Plot data

# MPS
MPS = []
for simu in DATA:
    rate = simu['interco']['rate']
    strength = simu['interco']['strength']
    MPS.append((rate, strength, simu['MPS']))
MPS.sort()

# plot MPS whole
X = linspace(0, 1, 10)
Y = linspace(0, 1, 10)
MESH = zeros((len(X), len(Y)))
for indx in xrange(len(X)):
    for indy in xrange(len(Y)):
        MESH[indx][indy] = MPS[indy*len(X) + indx][2]['whole']

CS = contourf(X, Y, MESH)
colorbar(CS)  # TODO mauvaise échelle...
show()



"""
DEPHASAGE

"""

PHI_ATTRS = (('paramset', '_v_attrs', 'Common', 'inter_conn_rate', 0, 1),
             ('paramset', '_v_attrs', 'Common', 'inter_conn_strength', 0, 1),
             ('results', '_v_attrs', 'phase_angles', 0, 1))
PHI = get_all_attrs(db, PHI_ATTRS)
PHI.sort()
X_PHI = list(set([i[0] for i in PHI]))
X_PHI.sort()
Y_PHI = list(set([i[1] for i in PHI]))
Y_PHI.sort()

Z_PHI = np.zeros((len(X_PHI), len(Y_PHI)))
for ind_rate in xrange(len(X_PHI)):
    for ind_strength in xrange(len(Y_PHI)):
        Z_PHI[ind_rate][ind_strength] = PHI[to1d(ind_rate, ind_strength, len(X_PHI))][2]

PHI_CS = contourf(X_PHI, Y_PHI, Z_PHI)
colorbar(PHI_CS)
show()

DB.close()
