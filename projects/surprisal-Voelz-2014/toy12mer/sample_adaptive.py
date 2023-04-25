# /usr/bin/env python

import random
import string
import math
import os, sys, copy, pickle

import scipy.sparse
import scipy.linalg
import scipy

from msmbuilder.MSMLib import *
from scipy.io import mmread, mmwrite

sys.path.append( os.path.join(os.environ['HPSANDBOXHOME'],'scripts') )
from HelperTools import *

sys.path.append('/Users/vince/scripts/msmbuilder/sandbox/vvoelz/Adaptive')
from ssaCalculatorVAV import ssaCalculator



usage = """Usage: sample_adaptive.py infile.mtx microstates.dat NumSamplesUniform NumSamplesAdaptive NumTrials outName

    This script will first perform an initial round of uniform sampling, 
    using NumSamplesUniform samples from each state (and the resampling disconnected states until ergodicity is reached)

    Then adaptive sampling based on Nina's eigenvalue sensitivity algorithm is performed, for NumTrials,
    sampling NumSamplesAdaptive from the most  

    Will write a series of files:
        tCounts.outName.NumSamplesAdaptive.0.mtx,
        ....
        tCounts.outName.NumSamplesAdaptive.(NumTrials-1).mtx

    

    Try: python sample_adaptive.py microT_AAAAAA.mtx microstates.dat 10 10000 10 adaptive
    print header

"""

if len(sys.argv) < 6:
    print usage
    sys.exit(1)
             
microTFn = sys.argv[1]
microstatesFn = sys.argv[2]
nsamples_uniform = int(sys.argv[3])
nsamples_adaptive = int(sys.argv[4])
ntrials = int(sys.argv[5])
outName = sys.argv[6]
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
if not os.path.exists('macrostates.dat'):
    fout = open('macrostates.dat','w')
    for i in range(len(mInfo.uniqueContactStates)):
        fout.write('%d\t%r\n'%(i,mInfo.uniqueContactStates[i]))
    fout.close()

############################
# generate increasing numbers of transition counts

NumStates = T.shape[0]
NumSamplesPerState = [nsamples_uniform for i in range(NumStates)]
NumTrials = ntrials

T = T.tolil()

jumpFn = microTFn+'.jumps'
if os.path.exists(jumpFn):
    print 'Reading jump probabilities from:', jumpFn, '...'
    fin = open(jumpFn,'r')
    jumps = pickle.load(fin)
    fin.close()
    print '...Done.'
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
TotalCounts = 0

NumContactStates = len(mInfo.uniqueContactStates)
C_macro  = scipy.sparse.lil_matrix((int(NumContactStates),int(NumContactStates)))
TotalMacroCounts = 0

# For the first step, do some uniform sampling of microstates until we get an ergodic network

ErgodicMacro = False  # is the macrostate transition matrix ergodic yet? 
while not ErgodicMacro:

    

    for i in np.array(NumSamplesPerState).nonzero()[0]:
        if i%1000 == 0:
            print 'Sampling', NumSamplesPerState[i], 'transitions from state', i, 'of', NumStates, '...'
        # look for possible transitions
        possible_j = jumps[i][0]
        jprobs = jumps[i][1]
        # sample counts from state i
        picks = draw_index(jprobs, n_picks=NumSamplesPerState[i])
        picked_counts = np.bincount(picks, minlength=len(possible_j))
        possible_macro_j = [mInfo.ContactStateIndices[j] for j in possible_j]
        for k in range(picked_counts.shape[0]):
            C[i,possible_j[k]] += picked_counts[k] 
            C_macro[mInfo.ContactStateIndices[i], possible_macro_j[k]] += picked_counts[k]
    TotalCounts += sum(NumSamplesPerState)
    TotalMacroCounts += sum(NumSamplesPerState)
    print 'Sampled', NumSamplesPerState[0], 'from each of 15037 states.  TotalCounts =', TotalCounts

    C_trimmed, Mapping = ErgodicTrim(scipy.sparse.csr_matrix(C_macro) )
    print 'C_trimmed.shape', C_trimmed.shape, 'C_macro.shape', C_macro.shape

    # Find all the non-ergodic states and sample only these in the next round 
    I_nonergodic = (Mapping == -1)
    MacroNumSamplesPerState = np.where(I_nonergodic, nsamples_uniform, 0) 
    NumSamplesPerState = [MacroNumSamplesPerState[mInfo.ContactStateIndices[i]] for i in range(NumStates)]

    if C_trimmed.shape == C_macro.shape:
        ErgodicMacro = True   


print 'C_macro', C_macro

# Next, use Nina's adaptive sampling algorithm to determine which macrostate should be sampled at each round 
for trial in range(NumTrials):

    lagTime = 1.
    s = ssaCalculator(lagTime, C_macro.todense(), evalList=[1,2], recommendationScheme = 'Nina', nNewSamples=nsamples_adaptive)
    EigInd = 0 
    MacroNumSamplesPerState = np.array(s.resamplingDistribution(EigInd))  # Recommended Samples
    print 'Adaptive sampling macrostate indices', MacroNumSamplesPerState.nonzero()[0]

    NumSamplesPerState = [MacroNumSamplesPerState[mInfo.ContactStateIndices[i]] for i in range(NumStates)]

    # sample microstates corresonding to the recommended macrostate
    for i in np.array(NumSamplesPerState).nonzero()[0]:
        if i%1000 == 0:
            print 'Sampling', NumSamplesPerState[i], 'transitions from state', i, 'of', NumStates, '...'
        # look for possible transitions
        possible_j = jumps[i][0]
        jprobs = jumps[i][1]
        # sample counts from state i
        picks = draw_index(jprobs, n_picks=NumSamplesPerState[i])
        picked_counts = np.bincount(picks, minlength=len(possible_j))
        possible_macro_j = [mInfo.ContactStateIndices[j] for j in possible_j]
        for k in range(picked_counts.shape[0]):
            C[i,possible_j[k]] += picked_counts[k]
            C_macro[mInfo.ContactStateIndices[i], possible_macro_j[k]] += picked_counts[k]
    TotalCounts += sum(NumSamplesPerState)
    TotalMacroCounts += sum(NumSamplesPerState)

    outfile = 'tCounts.micro.%s.%d.%d.mtx'%(outName,nsamples_adaptive,trial)
    print TotalCounts, 'samples. Writing', outfile
    mmwrite(outfile, C)

    outfile = 'tCounts.macro.%s.%d.%d.mtx'%(outName,nsamples_adaptive,trial)
    print TotalCounts, 'samples. Writing', outfile
    mmwrite(outfile, C_macro)






