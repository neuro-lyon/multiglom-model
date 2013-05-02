# -*- coding:utf-8 -*-
import matplotlib.pyplot as py
from scipy import *

# Path for parameter files must be created with "__init__.py" in it
paramfile_path='runs/one_glom/'
nproc=12
db_filename="db_one_glom_N100_sig003_highI.h5"


# General controls
make_sim=True
plot_results=True

# Fine simulation controls
clean_dir=make_sim
gen_param=make_sim
run_simul=make_sim
collect_db=make_sim

if clean_dir:
    import os
    list_files=os.listdir(paramfile_path)
    for f in list_files:
        if f!="__init__.py":
            os.remove(paramfile_path+"/"+f)
    print "Working directory cleaned"

if gen_param:
    # First generate parameter files

    from utils import *

    d = {('Common', 'N_subpop'): 
            {'range': [1],
             'unit': 1},
            ('Common', 'N_mitral'): 
            {'range': [100],
             'unit': 1},
             ('Common', 'simu_length'): 
            {'range': [2000],
             'unit': msecond},
            ('Input', 'g_Ein0'):
            {'range': linspace(start=0.1, stop=7., num=40),
             'unit': siemens*meter**-2},
            ('Input', 'sigma_Ein'):
            {'range': [0.03],
             'unit': siemens*meter**-2},
            ('Synapse', 'g_I'):
            {'range': [25.],
             'unit': siemens*meter**-2},
    }

    gen_parameters('paramsets/std_beta_1glom.py', d, paramfile_path)
    print "Params generated"

if run_simul:
    # Second run simulations in parallel


    from run_multip import *

    # Get all the parameter set files
    psfiles = get_pset_files(paramfile_path)

    # Run the simulation with each parameter set
    pool = Pool(processes=nproc)
    args = []
    for ind, psfile in enumerate(psfiles):
        args.append((psfile, ind, len(psfiles)))
        #~ new_simu(args[-1])

    print "Start simul"
    pool.map(new_simu, args)
    print "End simul"
    
if collect_db:
    # Third group result files in one hdf5 file
    from h5manager import *

    init_data_h5(db_filename)
    collect_h5_to_db(paramfile_path, db_filename, output=True)
    print "Data collected"
    
    
if plot_results:    
    # Fourth plot results
    from h5manager import *
    from plotting import *
    DB = tables.openFile(db_filename)

    # Plot full results for a list of simul
    list_plot=[0,0.1,0.2] # fraction of g_Ein0_max for plotted simul
    ATTRS = (('paramset', '_v_attrs', 'Input', 'g_Ein0'),
                 ('results', 'spikes_it'),
                 ('results', 's_granule'),
                 ('results', 's_syn_self'))
    res=get_all_attrs(DB, ATTRS)
    all_gEin0=array([r1 for r1,r2,r3,r4 in res])
    dt=DB.listNodes("/")[0].paramset._v_attrs['Common']['simu_dt']
    for frac in list_plot:
        # Find index of simul to plot
        g_Ein0_max=all_gEin0.max()
        ind=(abs(all_gEin0-frac*g_Ein0_max)).argmin()
        # Plot raster
        spit=res[ind][1]
        py.figure()
        py.scatter(spit[1,:],spit[0,:])
        py.title("g_Ein0 = "+str(all_gEin0[ind]))
        
        # Plot granule signals
        gr_s=res[ind][2][:]
        gr_s_syn_self=res[ind][3][:]
        times=arange(gr_s.size)*dt
        granule_pop_figure(gr_s, gr_s_syn_self, times, dt)


    # Plot synthtic results
    # Collect, g_Ein0, firing rate, network osc freq
    ATTRS = (('paramset', '_v_attrs', 'Input', 'g_Ein0'),
                 ('results', 'spikes_it'),
                 ('results', '_v_attrs', 'FFTMAX'))
    res=get_all_attrs(DB, ATTRS)

    simu_length=DB.listNodes("/")[0].paramset._v_attrs['Common']['simu_length']
    N_mitral=DB.listNodes("/")[0].paramset._v_attrs['Common']['N_mitral']
    start_time=simu_length/2.

    all_gEin0=[r1 for r1,r2,r3 in res]
    all_rates=[(r2[1,:]>=start_time).sum()/N_mitral/(simu_length-start_time) for r1,r2,r3 in res]
    all_freqs=[r3['mean'] for r1,r2,r3 in res]

    #~ print all_gEin0
    #~ print all_rates
    #~ print all_freqs

    fig=py.figure()
    ax=fig.add_subplot(121)
    ax.plot(all_gEin0,all_rates,'o',label="Firing rates")
    ax.plot(all_gEin0,all_freqs,'s',label="Oscillation freq")
    ax.set_xlabel("Input g_Ein0")
    ax.set_ylabel("Rate or freq (Hz)")
    ax.legend()

    ax=fig.add_subplot(122)
    ax.plot(all_rates,all_freqs,'o')
    ax.set_xlabel("Avg firing rate")
    ax.set_ylabel("Network freq")
    xl=ax.get_xlim()
    yl=ax.get_ylim()
    coord=min(xl[1],yl[1])
    coord2=min(xl[1],2*yl[1])
    ax.plot([0,coord],[0,coord],'--')
    ax.plot([0,coord2],[0,2*coord],'--')

    DB.close()
    py.show()

