#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from brian import *
from pylab import figure, plot, hist, show

def brun():
    eqs="""
    dZ/dt =  xi*second**-.5 : 1
    """

    nn   = NeuronGroup(1, model=eqs)
    print nn._eqs
    
    mon = {}
    for par in ['Z']:
        mon[par] = StateMonitor(nn, par, record=True)
    
    defaultclock.dt=0.01*msecond
    print 'dt =', defaultclock.dt
    
    run(500*msecond,report='text')
    
    for par in mon:
        figure()
        mon[par].plot()
        ylabel(par)
    
    print "Z SD noise :",mon['Z'][0].std()
    res = mon['Z'][0].std()

    show()

    defaultclock.reinit()
    clear(erase=True, all=True)

    return res

brun()
