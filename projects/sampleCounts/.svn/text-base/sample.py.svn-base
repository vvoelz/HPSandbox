#! /usr/bin/env python

import random
import string
import math
import os, sys, copy, pickle

import scipy.sparse
import scipy.linalg
import scipy

from msmbuilder.MSMLib import *

from scipy.io import mmread, mmwrite

sys.path.append('../')
from HelperTools import *

sys.path.append('/Users/vince/scripts/msmbuilder/sandbox/vvoelz')
from ssaCalculatorVAV import ssaCalculator


usage = """Usage: sample.py infile.mtx microstates.dat"""

if len(sys.argv) < 3:
    print usage
    sys.exit(1)
             
microTFn = sys.argv[1]
microstatesFn = sys.argv[2]
Adaptive = True

################
# Main program

# Read in transition matrix
if os.path.exists(microTFn):
    T = mmread(microTFn)
    print T
else:
    print "Can't find file:", microTFn, '...exiting'
    sys.exit(1)


# Read in microstate info
mInfo = MicrostateInfo(microstatesFn)

print 'mInfo', mInfo

############################
# generate increasing numbers of transition counts

NumStates = T.shape[0]
NumSamplesPerState = [1 for i in range(NumStates)]
NumTrials = 1000

T = T.tolil()

jumpFn = microTFn+'.jumps'
if os.path.exists(jumpFn):
    print 'Reading jump probabilities from:', jumpFn, '...'
    fin = open(jumpFn,'r')
    jumps = pickle.load(fin)
    fin.close()
else:
    print 'Making a dictionary of jump probabiilities...'
    jumps = {}   # {i, (possible_j, jprobs)}

    for i in range(NumStates):
        
        # look for possible transitions
        possible_j = T[i,:].nonzero()[1]
        #print i, '-->', 'possible_j', possible_j
        jprobs = np.array(T[i,possible_j].todense())
        jprobs = jprobs.reshape(jprobs.shape[1],)
        jumps[i] = (possible_j, jprobs)
        if i%100 == 0:
            print i, jumps[i]

    print 'pickling jumps to', jumpFn, '...'
    fout = open(jumpFn,'w')
    pickle.dump(jumps,fout)
    fout.close()
    print 'Done.'
    
    
C = scipy.sparse.lil_matrix((int(NumStates),int(NumStates))) 

#print '#ncounts\tstability\tslowest\tnext-slowest'
print '#ncounts\tstability\tslowest(tau1)\tnext-slowest(tau2)\tvar(eval1)\tresampling(e1)\tvar(eval2)\tresampling(e2)'

TotalCounts = 0
for trial in range(NumTrials):

    for i in range(NumStates):

        if NumSamplesPerState[i] > 0:

          if i%100 == 0: print 'resampling state', i 
        
          # look for possible transitions
          possible_j = jumps[i][0]
          jprobs = jumps[i][1]
        
          # sample counts from state i
          for k in range(NumSamplesPerState[i]):
              j = possible_j[draw_index(jprobs)]
              C[i,j] += 1
            
    TotalCounts += sum(NumSamplesPerState)
    
    T = EstimateTransitionMatrix(C)
    if (1): 
        # get the Implied Timescales
        EigAns, result = getTimescalesFromTransitionMatrix(T)
        EnoughCounts = True
    else:
        EnoughCounts = False
    
    if EnoughCounts: 
        # calculate the stability of the native state
        INative = np.array(mInfo.ContactStateIndices)==mInfo.NativeContactStateIndex
        stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )

        # print stability and timescale info 
        print '%d\t%8.3f\t%16.8f\t%16.8f'%(TotalCounts, stability, result[0,1], result[1,1])
    else:
        print '%d\tnot enough counts yet...'%(TotalCounts)


    # revise our estimates of how many samples to take from each state
    if Adaptive:
        lagTime = 1
        TotalSamples = sum(NumSamplesPerState) # this is the number of new samples we want to make each time'
        # print 'TotalSamples', TotalSamples
        s = ssaCalculator(lagTime, C, nNewSamples=TotalSamples, evalList=[1,2], recommendationScheme = 'Nina')
        #print 's.qlist', s.qlist
        #print 's.varianceContributions', s.varianceContributions
        for EigInd in range(len(s.evalList)):
            print s.uncertainty_variance(EigInd),
            #print 's.resamplingDistribution', 'for eigenvalue index', s.evalList[EigInd]
            RecommendedSamplesPerState = list(s.resamplingDistribution(EigInd))
            print np.argmax(RecommendedSamplesPerState)
            if EigInd == 0:
                NumSamplesPerState = RecommendedSamplesPerState



sys.exit(1)


    

