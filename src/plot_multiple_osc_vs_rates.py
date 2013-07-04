# -*- coding:utf-8 -*-
import matplotlib.pyplot as py
from matplotlib.cm import get_cmap
from scipy import *
from h5manager import *

"""
Script to plot osc vs rate curve coming from different set of simulations
"""

list_files=[['db_one_glom_N100_sig035_gI20.h5','db_one_glom_N100_sig035_gI20_gE1_4.h5','db_one_glom_N100_sig035_gI20_gE3_5.h5']]
                
def get_data(db_filename):
    
    DB = tables.openFile(db_filename)
 
    # Collect, g_Ein0, firing rate, network osc freq
    ATTRS = (('paramset', '_v_attrs', 'Input', 'g_Ein0'),
                 ('results', 'spikes_it'),
                 ('results', '_v_attrs', 'FFTMAX'))
    res=get_all_attrs(DB, ATTRS)

    simu_length=DB.listNodes("/")[0].paramset._v_attrs['Common']['simu_length']
    N_mitral=DB.listNodes("/")[0].paramset._v_attrs['Common']['N_mitral']
    start_time=simu_length/2.

    all_gEin0=[r1 for r1,r2,r3 in res]
    all_rates=[1.*(r2[1,:]>=start_time).sum()/N_mitral/(simu_length-start_time) for r1,r2,r3 in res]
    all_freqs=[r3['mean'] for r1,r2,r3 in res]

    DB.close()

    return all_gEin0,all_rates,all_freqs
    
    
n_rows=len(list_files)
fig=py.figure()

for ind,sublist_files in enumerate(list_files):
    
    ax=fig.add_subplot(n_rows,2,2*ind+1)
    cmap = get_cmap('jet', len(sublist_files))
    
    for num,db_filename in enumerate(sublist_files):
        all_gEin0,all_rates,all_freqs=get_data(db_filename)
        
        ax.scatter(all_rates,all_freqs,color=cmap(num),label=db_filename)
    ax.set_xlabel("Avg firing rate")
    ax.set_ylabel("Network freq")
    xl=ax.get_xlim()
    yl=ax.get_ylim()
    coord=min(xl[1],yl[1])
    coord2=min(xl[1],2*yl[1])
    ax.plot([0,coord],[0,coord],'--')
    ax.plot([0,coord2],[0,2*coord],'--')
    ax.set_xlim(0,40)
    ax.set_ylim(0,100)
    
    ax.legend()

    ax=fig.add_subplot(n_rows,2,2*ind+2)
    
    for num,db_filename in enumerate(sublist_files):
        all_gEin0,all_rates,all_freqs=get_data(db_filename)
        
        ax.scatter(all_gEin0,all_rates,color=cmap(num),label=db_filename)
    ax.set_xlabel("g_Ein0")
    ax.set_ylabel("Network freq")
    ax.set_xlim(0,5.2)
    ax.set_ylim(0,100)
    
py.show()