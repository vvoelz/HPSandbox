#! /usr/bin/env python

import os, sys, copy, pickle, random, string, math

import scipy.sparse
import scipy.linalg
import scipy

from msmbuilder.MSMLib import *

from scipy.io import mmread, mmwrite

sys.path.append( os.path.join(os.environ['HPSANDBOXHOME'],'scripts') )
from HelperTools import *


usage = """Usage: sample_traj.py infile.mtx microstates.dat NumSamples NumTrials outName

    Will write a series of files:
        tCounts.outName.NumSamples.0.mtx,
        ....
        tCounts.outName.NumSamples.(NumTrials-1).mtx

    Try: python sample_traj.py microT_AAAAAA.mtx microstates.dat 10000 10 test0 
"""

if len(sys.argv) < 6:
    print usage
    sys.exit(1)
             
microTFn = sys.argv[1]
microstatesFn = sys.argv[2]
nsamples = int(sys.argv[3])
ntrials = int(sys.argv[4])
outName = sys.argv[5]
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
    
state = 0    
C = scipy.sparse.lil_matrix((int(NumStates),int(NumStates)))
TotalCounts = 0

NumContactStates = len(mInfo.uniqueContactStates)
C_macro  = scipy.sparse.lil_matrix((int(NumContactStates),int(NumContactStates)))

obsfile = open('obs.traj.numcontacts.%s.dat'%outName,'w')
header = '#step\tnumContacts\tmacroindex\tmacrostate'
obsfile.write(header+'\n')
print header

for trial in range(NumTrials):

    for i in range(nsamples):

        #if i%1000 == 0: print 'Sampling', i, 'of', nsamples, ': state', state
        
        # look for possible transitions
        possible_j = jumps[state][0]
        jprobs = jumps[state][1]
        
        # sample counts from state i
        newstate = possible_j[draw_index(jprobs, n_picks=1)[0]]
        C[state, newstate] += 1
        contactState, newContactState = mInfo.ContactStateIndices[state], mInfo.ContactStateIndices[newstate]
        C_macro[contactState, newContactState] += 1

        TotalCounts += 1

        # print 'stair plot' number-of-contacts trajectory
        if newContactState != contactState:
            contactList = mInfo.uniqueContactStates[contactState]
            obsline = '%d\t%d\t%d\t%r'%(TotalCounts, len(contactList), contactState, contactList)
            obsfile.write(obsline+'\n')
            print obsline

            newContactList = mInfo.uniqueContactStates[newContactState]
            obsline = '%d\t%d\t%d\t%r'%(TotalCounts, len(newContactList), newContactState, newContactList)
            obsfile.write(obsline+'\n')
            print obsline

        state = newstate

    # Write micro count matrix
    outfile = 'tCounts.traj.micro.%s.%d.%d.mtx'%(outName,nsamples,trial) 
    print TotalCounts, 'samples. Writing', outfile
    mmwrite(outfile, C)

    outfile = 'tCounts.traj.macro.%s.%d.%d.mtx'%(outName,nsamples,trial)
    print TotalCounts, 'samples. Writing', outfile
    mmwrite(outfile, C_macro)

del C
    

