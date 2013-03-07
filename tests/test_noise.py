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
I_input = g*V: amp*meter**(-2)
"""

# NeuronGroup initialisation
nn   = NeuronGroup(1, model=eqs)
nn.g = g0
nn.V = 1

# Monitoring
defaultclock.dt = 0.01*ms
print 'NeuronGroup dt=', nn.clock.dt,
monit = StateMonitor(nn, 'g', record=True)

# Running simulation
print 'defaultclock dt=', defaultclock.dt,
run(200*msecond)

print 'g noise SD =', monit[0].std()
