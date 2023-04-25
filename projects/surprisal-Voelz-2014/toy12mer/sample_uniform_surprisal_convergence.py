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

sys.path.append( os.path.join(os.environ['HPSANDBOXHOME'],'scripts') )
from HelperTools import *

sys.path.insert(0, '/Users/vince/git/scripts/belman')
from Surprisal import *


usage = """Usage: sample.py infile.mtx NumSamplesPerState NumTrials outName

    Will write a series of files:
        tCounts.outName.NumSamplesPerState.0.mtx,
        ....
        tCounts.outName.NumSamplesPerState.(NumTrials-1).mtx

    Try: python sample_uniform_surprisal_convergence.py macroT_AAAAAA.mtx 10 100 surpconv
"""

if len(sys.argv) < 5:
    print usage
    sys.exit(1)
             
microTFn = sys.argv[1]
nsamples = int(sys.argv[2])
ntrials = int(sys.argv[3])
outName = sys.argv[4]
Adaptive = True


def uniform(jumps, nsamples, C_init=None):
    """Given a jumps dictionary and number of samples per state,
    do uniform sampling of transition counts and return a sparse count matrix."""

    NumStates = len(jumps.keys())
    if C_init != None:
        C = C_init.copy()
    else:
        C = scipy.sparse.lil_matrix((int(NumStates),int(NumStates)))

    for i in range(NumStates):
        print 'Uniform sampling for state', i
        # look for possible transitions
        possible_j = jumps[i][0]
        jprobs = jumps[i][1]

        # sample counts from state i
        for k in range(NumSamplesPerState[i]):
            j = possible_j[draw_index(jprobs)][0]
            C[i,j] += 1

    return C

def uniform_ergodic(jumps, nsamples, C_init=None):
    """Given a jumps dictionary and number of samples per state,
    do uniform sampling of transition counts and return a sparse count matrix.

    If nsamples of uniform sampling is *not* ergodic, keep sampling it until it is"""

    NumStates = len(jumps.keys())
    C = scipy.sparse.lil_matrix((int(NumStates),int(NumStates)))

    for i in range(NumStates):
        print 'Uniform sampling for state', i
        # look for possible transitions
        possible_j = jumps[i][0]
        jprobs = jumps[i][1]

        # sample counts from state i
        for k in range(NumSamplesPerState[i]):
            j = possible_j[draw_index(jprobs)][0]
            C[i,j] += 1

        while C[i,i] == C[i,:].sum():
            j = possible_j[draw_index(jprobs)][0]
            C[i,j] += 1

    return C





################
# Main program

# Read in transition matrix
if os.path.exists(microTFn):
    T = mmread(microTFn)
    print T
else:
    print "Can't find file:", microTFn, '...exiting'
    sys.exit(1)



############################
# generate increasing numbers of transition counts

NumStates = T.shape[0]
NumSamplesPerState = [nsamples for i in range(NumStates)]
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


# get a reference count matrix
Cref = uniform_ergodic(jumps, nsamples)

D_samples = []

for trial in range(NumTrials):

    C = uniform(jumps, nsamples, C_init=Cref) 
    print 'trial', trial, C

    # get estimated equilibrium populations
    T = EstimateTransitionMatrix(C)
    EigAns = GetEigenvectors(T,5)
    pi = EigAns[1][:,0]
    pi = pi/pi.sum() # normalize
    print 'T', T
    print 'pi', pi

    print 'C', C
    print 'Cref', Cref

    C = C.tocsr()
    Cref = Cref.tocsr() 

    print Cref.row

    obj = SurprisalCalculator(C, Cref)
    D_i = []  # a list of D_i for each state
    for i in range(NumStates):
        n1, n2 = obj.prepare_count_arrays(i)
        s_i = obj.calculate_surprisal(n1, n2, weighted=True)
        D.append(pi[i]*s_i)
    D = D_i.sum()

    D_samples.append(D)

    TotalCounts = C.sum()
    print 'Sampled', trial*nsamples, 'from each of', NumStates, 'states.  TotalCounts =', TotalCounts

print 'D_samples', D_samples
