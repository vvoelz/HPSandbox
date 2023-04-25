#! /usr/bin/env python

import random
import string
import math
import os, sys, copy, pickle, glob

import numpy as np

import scipy.sparse
import scipy.linalg
import scipy

from scipy.io import mmread, mmwrite

from msmbuilder.MSMLib import *

from scipy.io import mmread, mmwrite

sys.path.append('../')
from HelperTools import *

def RemoveEmptyRows(C):
    # removes empty rows and cols from a COO sparse matrix
    rows = C.row
    cols = C.col
    data = C.data

    uniq = list(np.unique(np.array(list(np.unique(rows)) + list(np.unique(cols)))))
    newrow = [uniq.index(r) for r in rows]
    newcol = [uniq.index(c) for c in cols]

    return scipy.sparse.coo_matrix((data,(newrow,newcol)), shape=(len(uniq),len(uniq)))




print '#ncounts\tstatesAfterTrim\tstability\tslowest(tau1)\tnext-slowest(tau2)'

#count_matrix_Fns = glob.glob("tCounts.traj.microT_AAAAAA.1500000.*.mtx")
#count_matrix_Fns = ["tCounts.traj.microT_AAAAAA.15000.%d.mtx"%i for i in range(100)]
#count_matrix_Fns = ["tCounts.uniform.microT_AAAAAA.1.%d.mtx"%i for i in range(1,100)]
count_matrix_Fns = ["tCounts.uniform.microT_AAAAAA.100.%d.mtx"%i for i in range(100)]

# count_matrix_Fns = ["tCounts.equil.microT_AAAAAA.15000.%d.mtx"%i for i in range(10,100)]

VERBOSE = False

CountType = 'uniform'
if count_matrix_Fns[0].count('traj') > 0:
    CountType = 'traj'
elif count_matrix_Fns[0].count('equil') > 0:
    CountType = 'equil'

for count_matrix_Fn in count_matrix_Fns:

    TotalCounts = int(count_matrix_Fn.split('.')[-2])*int(count_matrix_Fn.split('.')[-3])
    if CountType == 'uniform':
        TotalCounts *= 15037

    D = mmread(count_matrix_Fn)

    if VERBOSE: print 'Doing ergodic trim...'
    C, MAP = ErgodicTrim(D)
    if VERBOSE: print 'Trimmed from', D.shape[1], 'to', C.shape[0], 'indices'


    T = EstimateTransitionMatrix(C)

    # get the Implied Timescales
    EigAns, result = getTimescalesFromTransitionMatrix(T, NumImpliedTimes= min(6,T.shape[0]-1))
    
    # calculate the stability of the native state
    #INative = np.array(mInfo.ContactStateIndices)==mInfo.NativeContactStateIndex
    INative = 11811
    if VERBOSE: print 'new index for native state:', INative, '-->', MAP[INative]
    #stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )
    #if TrajCounts:
    #    stability = D.sum(1)[INative]/D.sum(1).sum(0)
    #else:
    if MAP[INative] != -1:
        stability = np.real( EigAns[1][MAP[INative],0]/EigAns[1][:,0].sum() )
    else:
        stability = 0.

    # print stability and timescale info 
    print '%d\t%d\t%8.3f\t%16.8f\t%16.8f'%(TotalCounts, T.shape[0], stability, result[0,1], result[1,1]), count_matrix_Fn

sys.exit(1)


