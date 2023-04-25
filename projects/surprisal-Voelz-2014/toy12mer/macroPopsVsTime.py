#!/usr/bin/env python

import os, sys
import numpy as np

from msmbuilder import MSMLib
from scipy.io import mmread
from scipy import loadtxt, savetxt

sys.path.append('../../scripts')
from HelperTools  import *

import matplotlib.pyplot as plt

usage = """Usage:  python macroPopsVsTime.py microT_AAAAAA.mtx macroT_AAAAAA.mtx microstates.dat   """

if len(sys.argv) < 3:
    print usage
    sys.exit(1)   

tProbsFn = sys.argv[1]
tProbsFn_macro = sys.argv[2]
microstatesFn = sys.argv[3]

T_micro = mmread(tProbsFn).tocsr()
print 'T_micro.shape', T_micro.shape

T_macro = mmread(tProbsFn_macro).todense() 
print 'T_macro.shape', T_macro.shape

m = MicrostateInfo(microstatesFn)
nstates = len(m.uniqueContactStates)

# Calculate time evolution from T_macro 
p = np.zeros( (nstates,) )
p[0] = 1.0
print 'p.shape', p.shape

nsteps = 50000
result = np.zeros( (nsteps, nstates) )
for step in range(nsteps):
    result[step,:] = p
    p =  (np.dot(p, T_macro)).flatten()

savetxt('macroPopsVsTime_fromTmacro.dat', result)

#####
print 'Calculate time evolution from T_micro ...'
nstates_micro = T_micro.shape[0]
T_microT = T_micro.transpose()
q = np.zeros( (nstates_micro,) )
q[0] = 1.0
print 'q.shape', q.shape

nsteps_micro = 50000
result_micro = np.zeros( (nsteps_micro, nstates) )

for step in range(nsteps_micro):
    for i in range(nstates):
        qselect = np.where(np.array(m.ContactStateIndices) == i, 1, 0)
        #print 'i', i, 'qselect', qselect
        p[0,i] = np.dot(q,qselect) 
    print 'step', step, 'p', p[0,0:5], '...'
    result_micro[step,:] = p.flatten()
    q =  T_microT.dot(q)

savetxt('macroPopsVsTime_fromTmicro.dat', result_micro)


fig = plt.figure()

ax = fig.add_subplot(1,2,1)
ax.plot( np.arange(nsteps), result)
ax.axis( [1, nsteps, 0.03, 1.0] )
ax.set_xscale('log')
ax.set_yscale('log')

ax = fig.add_subplot(1,2,2)
ax.plot( np.arange(nsteps_micro), result_micro)
ax.axis( [1, nsteps, 0.03, 1.0] )
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()



