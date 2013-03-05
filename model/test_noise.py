#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from brian import *

# Parameters
g0     = 0.5*siemens*meter**-2.
sigma  = 0.002*siemens*meter**(-2.)*second**(0.5)
tg     = 3*msecond
g_Ein  = 1.*siemens*meter**(-2)

# Model equations
eqs="""
dV/dt   = 0*volt/second: volt
dg/dt   = (-g + g0 + sigma * xi)/tg : siemens*meter**(-2)
I_input = g*volt: amp*meter**(-2)
"""

# NeuronGroup initialisation
nn   = NeuronGroup(1, model=eqs)
nn.g = g0
nn.V = 1

# Monitoring
monit = StateMonitor(nn, 'g', record=True)

# Running simulation
defaultclock.dt = 0.001*msecond
print 'dt =', defaultclock.dt

run(1000*msecond, report='text')

print "g SD noise :", monit[0].std()
