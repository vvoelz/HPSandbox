#!/usr/bin/env python

import os, sys
import numpy as np

from msmbuilder import MSMLib
from scipy.io import mmread
from scipy import loadtxt, savetxt

sys.path.append('../../scripts')
from HelperTools  import *

import matplotlib.pyplot as plt

print 'Loading macroPopsVsTime_fromTmacro.dat ...'
pops_macro = loadtxt('macroPopsVsTime_fromTmacro.dat')

print 'Loading macroPopsVsTime_fromTmicro.dat ...'
pops_micro = loadtxt('macroPopsVsTime_fromTmicro.dat')

m = MicrostateInfo('microstates.dat')
#print m.uniqueContactStates

# print top populations in order
Isort = np.argsort(-pops_macro[-1,:])
ntop = 10
for i in range(ntop):
   print pops_macro[-1,Isort[i]], Isort[i], m.uniqueContactStates[Isort[i]]
    
