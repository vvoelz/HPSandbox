#! /usr/bin/env python

import random
import string
import math
import os, sys, copy, pickle 
import numpy as np

import scipy.sparse
import scipy.linalg
import scipy

from scipy.io import mmread, mmwrite

from msmbuilder.MSMLib import *

from scipy.io import mmread, mmwrite

sys.path.append('../')
from HelperTools import *

sys.path.append('/Users/vince/scripts/msmbuilder/sandbox/vvoelz')
from ssaCalculatorVAV import ssaCalculator

toy12mer = True

if toy12mer:

    microTFn = '../build12merTmat/microT_AAAAAA.mtx'
    jumpFn = '../build12merTmat/microT_AAAAAA.mtx.jumps'

    # Read in the transition matrix
    T = mmread(microTFn)
    print 'T', T

    print 'Reading jump probabilities from:', jumpFn, '...'
    fin = open(jumpFn,'r')
    jumps = pickle.load(fin)
    fin.close()
    print '...Done.'

    NumTrials = 1000
    NumStates = T.shape[0]   
    NumSamplesPerState = [10 for i in range(NumStates)]
    Adaptive = True


else:
    NumTrials = 1000
    NumStates = 10   
    NumSamplesPerState = [10 for i in range(NumStates)]
    Adaptive = True

    # make a random transition matrix
    TFn = 'rand%d.mtx'%NumStates
    if os.path.exists(TFn):
        T = mmread(TFn)
    else:
        T = EstimateTransitionMatrix( np.random.random( (NumStates,NumStates) ) )
        mmwrite(TFn, T, precision=32)
    print 'T', T

    # make an empty count matrix
    C = scipy.sparse.lil_matrix( (NumStates,NumStates) ) 

    # build jump probabilities from the transition matrix
    jumps = {}
    for i in range(NumStates):
        jumps[i] = (T[i,:].nonzero()[0], T[i,:])
    print 'jumps', jumps

# make an empty count matrix
C = scipy.sparse.lil_matrix( (NumStates,NumStates) )


print '#ncounts\tstability\tslowest(tau1)\tnext-slowest(tau2)\tvar(eval1)\tresample?(e1)\tvar(eval2)\tresample?(e2)'
TotalCounts = 0
for trial in range(NumTrials):

    #print 'Trial', trial, '- sampling', NumSamplesPerState[0], 'for each of', NumStates, 'states.'
    for i in range(NumStates):

        if i%1000 == 0: print 'sampling state', i, 'of', NumStates, '...'
        
        # look for possible transitions
        possible_j = jumps[i][0]
        jprobs = jumps[i][1]
        
        # sample counts from state i
        for k in range(NumSamplesPerState[i]):
            j = possible_j[draw_index(jprobs)]
            C[i,j] += 1

        # symmetrize the count matrix!  If not, we'll get a pair of a+bi, a-bi eigenvalues, so that the real components
        # yield two identical eigenvectors.
        #for j in range(i): 
        #    C[i,j] = (C[i,j]+C[j,i])/2.0
   
            
    TotalCounts += sum(NumSamplesPerState)
    
    T = EstimateTransitionMatrix(C)
    if (1): 
        # get the Implied Timescales
        EigAns, result = getTimescalesFromTransitionMatrix(T, NumImpliedTimes = 6)
        EnoughCounts = True
    else:
        EnoughCounts = False
    #print 'EnoughCounts?', EnoughCounts
    
    if EnoughCounts: 
        # calculate the stability of the native state
        #INative = np.array(mInfo.ContactStateIndices)==mInfo.NativeContactStateIndex
        #stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )
        stability = 0.0  # for testing

        # print stability and timescale info 
        print '%d\t%8.3f\t%16.8f\t%16.8f'%(TotalCounts, stability, result[0,1], result[1,1]),
    else:
        print '%d\tnot enough counts yet...'%(TotalCounts),


    # revise our estimates of how many samples to take from each state
    if Adaptive:
        lagTime = 1
        TotalSamples = sum(NumSamplesPerState) # this is the number of new samples we want to make each time'
        #print 'TotalSamples', TotalSamples
        s = ssaCalculator(lagTime, C.todense(), nNewSamples=TotalSamples, evalList=[1,2], recommendationScheme = 'Nina')
        #print 's.varianceContributions', s.varianceContributions
        for EigInd in range(len(s.evalList)):
            #print s.uncertainty_variance(EigInd),
            if EigInd == 0:
                RecommendedSamplesPerState = list(s.resamplingDistribution(EigInd))
                NumSamplesPerState = RecommendedSamplesPerState
        if (0):
            print 
            print 's.weights', np.array(s.weights)
            print 's.qlists[EigInd=0]', s.qlists[0]
            print 's.resamplingDistribution', 'for eigenvalue index', s.evalList[0], repr(RecommendedSamplesPerState).replace(' ','').strip()

    print '----------'

sys.exit(1)


