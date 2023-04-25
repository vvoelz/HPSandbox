#!/usr/bin/env python

import os, sys
import numpy as np

from msmbuilder import MSMLib
from scipy.io import mmread
from scipy import loadtxt, savetxt

usage = """Usage:  python getPops.py tProbsFn.mtx
    will write tProbsFn.Populations.dat"""

if len(sys.argv) < 2:
    print usage
    sys.exit(1)   

tProbsFn = sys.argv[1]
TC = mmread(tProbsFn)
EigAns = MSMLib.GetEigenvectors(TC,5)
Populations = EigAns[1][:,0]

popsFn = os.path.basename(tProbsFn).replace('.mtx','.Populations.dat')
savetxt(popsFn, Populations)

