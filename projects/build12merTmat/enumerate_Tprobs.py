#! /usr/bin/python

usage = """
enumerate.py <configfile>

Try:  enumerate_Tprobs.py enumerate.conf

This program will read in an HP chain specified in the configure file,
and perform a full enumeration of conformational space.

The problem tablulates:

    1) the density of states (in energies/contacts)
    
    2) the number density of unique contact states, i.e. disjoint collections
       of microscopic conformations all sharing a unique set of interresidue contacts. 

These values are printed as output.

"""


import sys
sys.path.append('../../hpsandbox')

from Config import *
from Chain import *
from Monty import *
from Replica import *
from Trajectory import *

import random
import string
import math
import os, copy

import scipy.sparse
import scipy.linalg
import scipy

from msmbuilder.MSMLib import *


g = random.Random(randseed)


if len(sys.argv) < 2:
    print('Usage:  enumerate_Tprobs.py <configfile>')
    sys.exit(1)
 
 
VERBOSE = 1

def nonsym(vec):
        """Many of the conformations are related by rotations and reflections.
        We define a "non-symmetric" conformation to have the first direction '0'
        and the first turn be a '1' (right turn)
        
        nonsym() returns 1 if the vec list is non-symmetric, 0 otherwise
        """
        if len(vec) > 0:
            # walk along chain until you get to the first non '0' vec:
            i = 0
            while (vec[i] == 0) & (i < len(vec) - 1) :
                i = i + 1
            if vec[0] == 0:
                if (vec[i] == 1) | (vec[i] == 0):
                    return(1)

        return(0)


def vec2coords(thisvec):
    """Convert a list of chain vectors to a list of coordinates (duples)."""
    tmp = [(0,0)]
    x = 0
    y = 0
    for i in range(0,len(thisvec)):
        if thisvec[i] == 0:
            y = y + 1
        if thisvec[i] == 1:
             x = x + 1
        if thisvec[i] == 2:
            y = y - 1
        if thisvec[i] == 3:
            x = x - 1
        tmp.append((x,y))
    return tmp

def viable(thisvec):
    """Return True if the chain in viable (self-avoiding), False if not."""
    return bool(viability(vec2coords(thisvec)))

def viability(thesecoords):
    """Return 1 if the chain coordinates are self-avoiding, 0 if not."""
    viable = 1
    for c in thesecoords:
        if thesecoords.count(c) > 1:
            viable = 0
            break

    return viable




if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(1)
    
    configfile = sys.argv[1]
    config = Config( filename=configfile)
    if VERBOSE: config.print_config()
    
    # create a single Replica
    replicas = [ Replica(config,0) ]
    
    traj = Trajectory(replicas,config)    # a trajectory object to write out trajectories

    nconfs = 0
    contact_states = {}        # dictionary of {repr{contact state}: number of conformations}
    contacts = {}               # dictionary of {number of contacts: number of conformations}


    microstates = []            # a list of microstate vec lists, in order
    microContactStates = []     # a list of contact states, in order, for each microstate
    microindex = {}             # dictionary of {"011121321": microstate index} lookup table
    


    #################
    #
    # This is a useful subroutine for enumerating all conformations of an HP chain
    #
    # NOTE: in order for this to work correctly, the initial starting vector must be [0,0,0,....,0]
    # 
    done = 0
    while not(done):

        if len(replicas[0].chain.vec) == replicas[0].chain.n-1:    
            if replicas[0].chain.viable:        
                if replicas[0].chain.nonsym():
            
                    # tally the number of contacts
                    state = replicas[0].chain.contactstate()
                    ncontacts = len(state)
                    if contacts.has_key(ncontacts) == False:
                        contacts[ncontacts] = 1
                    else:
                        contacts[ncontacts] = contacts[ncontacts] + 1

                    # tally the contact state
                    this_state_repr = repr(state)
                    if contact_states.has_key(this_state_repr) == False:
                        contact_states[this_state_repr] = 1
                    else:
                        contact_states[this_state_repr] = contact_states[this_state_repr] + 1

                    # tally the number of conformations
                    nconfs = nconfs + 1

                    microstates.append(copy.copy(replicas[0].chain.vec))
                    microindex[ replicas[0].chain.vec2str(replicas[0].chain.vec) ] = len(microstates)-1 
                    microContactStates.append( state )

                    # write to trajectory
                    # if (nconfs % config.TRJEVERY) == 0:
                    #        traj.queue_trj(replicas[0])
                    # print progress
                    # if (nconfs % config.PRINTEVERY) == 0:
                    #    print('%-4d confs  %s'%(nconfs,replicas[0].chain.vec))
    
                done = replicas[0].chain.shift()
                    
            else:
                done = replicas[0].chain.shift()

        else:
            if replicas[0].chain.viable:
                replicas[0].chain.grow()
            else:
                done = replicas[0].chain.shift()

        if replicas[0].chain.vec[0] == 1:    # skip the other symmetries
            break    
    #
    #
    #################
        
    
    # write the last of the trj and ene buffers
    # and close all the open trajectory file handles
    traj.cleanup(replicas)
    
    # print out the density of contact states
    print()
    print('DENSITY of CONTACT STATES:')
    print('%-40s %s'%('contact state','number of conformations'))
    for state in contact_states.keys():
        print('%-40s %d'%(state, contact_states[state]))
    
    # print out the density of states (energies)
    print()
    print('DENSITY of STATES (in energies/contacts):')
    print('%-20s %-20s %s'%('number of contacts','energy (kT)','number of conformations'))
    for c in contacts.keys():
        print('%-20d %-20d %d'%(c,config.eps*c,contacts[c]))
    print()
    print('at T = %4.1f K'%config.T)

    # write a file describing each microstate
    microstateFn = 'microstates.dat'
    print('Writing %s ....'%microstateFn, end='') 
    fout = open(microstateFn, 'w')
    for i in range(len(microstates)):
        fout.write('%d\t%s\t%d\t%s\n'%(i, repr(microstates[i]), len(microContactStates[i]), repr(microContactStates[i])))
    fout.close()
    print('Done.')


    ############
    # For each microstate, compile a transition matrix for all possible moves

    from scipy.io import mmread, mmwrite
    microTFn = 'microT.mtx'
    if os.path.exists(microTFn):
         T = mmread(microTFn)
         print(T)
    else:

        NumStates = len(microstates)
        C=scipy.sparse.lil_matrix((int(NumStates),int(NumStates)))  # Lutz: why are we using float for count matrices?

        for i in range(len(microstates)):

            if i%100 == 0:
                print('computing (iso-energetic) transitions from state', i, 'of', NumStates)

            # iterate through all possible moves
            for vecindex in range(0,replicas[0].chain.n):
              for direction in [-1,1]:
                for case in [0,1,2]:

                  vec = microstates[i]
                  nextvec = copy.copy(microstates[i])

                  # if possible, 1/3 of the time do a three-bead flip (dirs must be different)
                  if (case == 0) and (vecindex < len(nextvec)-1):
                      tmp1 = nextvec[vecindex]
                      tmp2 = nextvec[vecindex+1]
                      if (tmp1 != tmp2):
                          nextvec[vecindex] = tmp2
                          nextvec[vecindex + 1] = tmp1
                      else:
                          ### default: do a rigid rotation
                          for v in range(vecindex,len(nextvec)):
                              nextvec[v] = (nextvec[v] + direction) % 4

                  # if possible, 1/3 of the time do a crankshft (1st and 3rd dirs must be different)
                  elif ((case == 2) and (vecindex < len(nextvec)-2)):
                      tmp1 = nextvec[vecindex]
                      tmp2 = nextvec[vecindex+2]
                      if (tmp1 != tmp2):
                          ### crankshaft move ###
                          nextvec[vecindex] = tmp2
                          nextvec[vecindex + 2] = tmp1
                      else:
                          ### default: do a rigid rotation
                          for v in range(vecindex,len(nextvec)):
                              nextvec[v] = (nextvec[v] + direction) % 4

                  # the remaining 1/3 of the time, do a rigid rotation
                  else:
                      ### default: do a rigid rotation
                      for v in range(vecindex,len(nextvec)):
                          nextvec[v] = (nextvec[v] + direction) % 4


                  ### Find microstate index we're transition to

                  # correct rotational symmetry
                  if nextvec[0] != 0:
                      for v in range(len(nextvec)):
                          nextvec[v] = (nextvec[v] - nextvec[0]) % 4

                  # correct mirror symmetry
                  if not nonsym(nextvec):
                      for v in range(len(nextvec)):
                          nextvec[v] = (-nextvec[v]) % 4

                  if viable(nextvec):
                      # add to the count matrix
                      j = microindex[replicas[0].chain.vec2str(nextvec)]
                      C[i,j] += 1.0

                      if (0): 
                          print('\t%s(%d)'%(replicas[0].chain.vec2str(vec), i),  '-->' , end='')
                          print('\t%s(%d)'%(replicas[0].chain.vec2str(nextvec), j))
                  else:
                      C[i,i] += 1.0

        # Calculate the eigen problem
        if (1):
            T = EstimateTransitionMatrix(C)
        else:
            weights = np.asarray(C.sum(axis=1)).flatten()
            D=scipy.sparse.dia_matrix((1./weights,0),C.shape).tocsr()
            T=D.dot(C)

        print('Writing microstate transition matrix to %s ...'%microTFn, end='')
        mmwrite(microTFn, T, precision=32)   # write the transition matrix to file (for next time)
        print('...Done.')


    print('Computing Implied Timescales....')
    LagTime = 1.   # one step
    NumImpliedTimes = 10
    EigAns=GetEigenvectors(T,NumImpliedTimes+1,Epsilon=1) #TJL: set Epsilon high, should not raise err here     

    # make sure to leave off equilibrium distribution
    lagTimes = LagTime*np.ones((NumImpliedTimes))
    impTimes = -lagTimes/np.log(EigAns[0][1:NumImpliedTimes+1])

    # save intermediate result in case of failure
    result = np.zeros((NumImpliedTimes, 2))
    result[:,0] = lagTimes
    result[:,1] = impTimes

    print(result)

