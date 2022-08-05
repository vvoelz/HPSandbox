#! /usr/bin/python

usage = """
mcrex.py <configfile>

Try:  mcrex.py mcrex.conf

This program will read in an HP chain and run parameters specified in the configure file,
and perform a replica exchange Monte Carlo simulation.

For the example "mcrex.conf", an 11-mer sequence is simulated, and the program ends when
the native conformation (contact state) is found.
A directory of results is output to directory ./mcrex_data
"""

import sys
sys.path.append('../')

from Config import *
from Chain import *
from Monty import *
from Replica import *
from Trajectory import *

import random
import string
import math
import os

g = random.Random(randseed)


if len(sys.argv) < 2:
    print usage
    sys.exit(1)
 
 
VERBOSE = True
    

    
##########################
# Functions 

def attemptswap(replicas,swaps,viableswaps): 

    if len(replicas) == 1:
       return 0 

    # Attempt a swap between replicas    
    if config.SWAPMETHOD == 'random pair':
        # pick pair at random
        r = g.random()
        i = min(int(r*config.NREPLICAS),(config.NREPLICAS-1))
        j = i
        while (j==i):
            s = g.random()
            j = min(int(s*config.NREPLICAS),(config.NREPLICAS-1))
        # make sure j > i (needed for the swap criterion below)
        if j < i:
            tmp = i
            i = j
            j = tmp
    
    elif config.SWAPMETHOD == 'neighbors':    
        # pick neighboring pair at random 
        r = g.random()
        i = min(int(r*(config.NREPLICAS-1)),(config.NREPLICAS-2))
        j = i+1
                
    else:
        print 'Swap method', config.SWAPMETHOD, 'unknown.  Exiting...'
        sys.exit(1)


    if (VERBOSE):
      print 'REX: Attempting swap between replicas',i,'and',j
     
    randnum = g.random()

    ### if proposing i-->j, 
    ### 
    ###   del = d(beta_j - beta_i) * d(E_j - E_i)
    ###
    ### (see Hansmann 1997)

    delfactor = ((1/(replicas[j].mc.k*replicas[j].mc.temp)) - (1/(replicas[j].mc.k*replicas[i].mc.temp)))
    delfactor = delfactor * (replicas[j].mc.energy(replicas[j].chain) - replicas[i].mc.energy(replicas[i].chain) )

    boltzfactor = math.exp(delfactor)
    if randnum < boltzfactor:

        # swap the ***temperatures***
        temp_i = replicas[i].mc.temp
        temp_j = replicas[j].mc.temp
        tempfromrep_i = replicas[i].mc.tempfromrep
        tempfromrep_j = replicas[j].mc.tempfromrep
        replicas[i].mc.temp = temp_j
        replicas[j].mc.temp = temp_i
        replicas[i].mc.tempfromrep = tempfromrep_j
        replicas[j].mc.tempfromrep = tempfromrep_i
        
        retval = 1
        
    else:
        retval = 0

    viableswaps[i] = viableswaps[i] + retval
    viableswaps[j] = viableswaps[j] + retval
    swaps[i] = swaps[i] + 1
    swaps[j] = swaps[j] + 1

    return retval







#####################################
# Main Program
    
if __name__ == '__main__':
    
    # load in config file
    configfile = sys.argv[1]
    config = Config()
    config.read_configfile(configfile)
    if VERBOSE:
        config.print_config()

    # look up the native contact state from a pre-compiled flat file     
    if config.STOPATNATIVE == 1:
        nativeclistfile = config.NATIVEDIR + '/' + config.HPSTRING + '.clist'
        fnative = open(nativeclistfile,'r')
        nativeclist = eval(fnative.readline())
        fnative.close()
        if VERBOSE:
            print 'NATIVE CLIST:',nativeclist
    
    
    # Make a list of config.NREPLICAS Replica objects
    replicas = []           # a list of Replica objects
    for i in range(0,config.NREPLICAS):
        replicas.append(Replica(config,i))
    
    # Each replica is just a container for the objects:
    #     replica[i].chain   <-- the HP chain
    #     replica[i].mc      <-- A Monte Carlo move set object that performs operations on the chain
    
    # a trajectory object to handle the work of writing trajectories
    traj = Trajectory(replicas,config)  
    if VERBOSE:
        print 'Trajectory REPFILES:',traj.repfiles_trj


    # Keep track of statistics for the replica exchange simulation, for each replica:

    ### Monte Carlo stats
    steps = []                 # number of total move steps attempted (viable or not)
    viablesteps = []           # number of viable steps
    acceptedsteps = []         # number of viable steps that were accepted under Metropolis criterion 
    moveacceptance = []        # fraction of 
    acceptance = []            # fraction of MC moves accepted

    ### Replica exchange stats 
    swaps = []                 # the number of attemped replica swaps
    viableswaps = []           # the number of viable swaps
    swap_acceptance = []       # the fraction of swaps that were accepted
    
    for i in range(0,config.NREPLICAS):
         steps.append(0)
         viablesteps.append(0)    
         acceptedsteps.append(0)    
         moveacceptance.append(0)    
         acceptance.append(0)    
         swaps.append(0)
         viableswaps.append(0)
         swap_acceptance.append(0)
    
    
    prodstep = 1
    foundnative = 0
    while (prodstep < (config.MCSTEPS+1))&(foundnative==0):
    
        # Run the replicas for a production cycle...
        for rep in range(0,config.NREPLICAS):

            ### Propose a new MC move
            if string.strip(config.MOVESET) == 'MS1':
                replicas[rep].mc.move1(replicas[rep])
            elif string.strip(config.MOVESET) == 'MS2':
                replicas[rep].mc.move2(replicas[rep])
            elif string.strip(config.MOVESET) == 'MS3':
                replicas[rep].mc.move3(replicas[rep])
            elif string.strip(config.MOVESET) == 'MS4':
                replicas[rep].mc.move4(replicas[rep])
            else: 
                print 'MC MOVESET=',config.MOVESET,'not supported!'
                print 'Exiting....'
                sys.exit(1)

            if replicas[rep].chain.nextviable == 1:     # count this move only if the chain is viable
                ### accept with metroplis criterion
                accept = replicas[rep].mc.metropolis(replicas[rep])
                # if accept is yes, the chain is updated (see Monty.py)
                if (accept):
                    acceptedsteps[rep] = acceptedsteps[rep] + 1
                # keep track of MC steps
                viablesteps[rep] = viablesteps[rep] + 1

            # keep track of total steps
            steps[rep] = steps[rep] + 1

            # write the trajectory with specified frequency
            if (steps[rep] % config.TRJEVERY) == 0:
                traj.queue_trj(replicas[rep])
            if (steps[rep] % config.ENEEVERY) == 0:
                traj.queue_ene(replicas[rep])

            # HAVE WE FOUND THE NATIVE YET? If so, stop
            if config.STOPATNATIVE == 1:    
                thisclist  = replicas[rep].chain.contactstate()
                if (nativeclist == thisclist):
                    foundnative = 1


        # After the production cycle,      
        if (prodstep % config.SWAPEVERY) == 0:

            ### ...after every production run, attempt a SWAP
            success = attemptswap(replicas,swaps,viableswaps)

            if (VERBOSE):
              if success:   
                outstring = str(prodstep)+'\tSwap successful!\treplica temps: ['
                for rep in range(0,config.NREPLICAS):
                  outstring = outstring + str(replicas[rep].mc.temp) + ' '
                print outstring +']'

              else:  
                print str(prodstep)+'\tUnsuccessful swap attempt.'


        # Print status
        if (prodstep % config.PRINTEVERY) == 0:

            # calc MC move acceptance
            for rep in range(0,len(steps)):
                if steps[rep] > 0:
                    moveacceptance[rep] = float(viablesteps[rep])/float(steps[rep])
                    if viablesteps[rep] > 0:
                        acceptance[rep] = float(acceptedsteps[rep])/float(viablesteps[rep])
                    else:
                        acceptance[rep] = 0.0
                else:
                    moveacceptance[rep] = 0.0
                    acceptance[rep] = 0.0

            # calc replica swap acceptance
            for rep in range(0,len(steps)):
                if swaps[rep] > 0:
                    swap_acceptance[rep] = float(viableswaps[rep])/float(swaps[rep])
                else:
                    swap_acceptance[rep] = 0.0

            # Output the status of the simulation
            print prodstep,'production steps'
            print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s '%('replica','temp(K)','viablesteps','steps','MCaccept','viableswaps','swaps','SWAPaccept')
            for rep in range(0,len(steps)):
                temp = replicas[rep].mc.temp
                print '%-12d %8.3f %-12d %-12d %-12s %-12d %-12d %-12s '%(rep,temp, viablesteps[rep],steps[rep],'%1.3f'%moveacceptance[rep],swaps[rep],viableswaps[rep],'%1.3f'%swap_acceptance[rep])
            if config.STOPATNATIVE == 1:
                print 'NATIVE CLIST:', nativeclist
            print '%-8s %-12s %-12s %-12s'%('replica','foundnative','contact state', 'chainvec')
            for replica in replicas:
                print '%-8d %-12d %-12s %-12s'%(replica.repnum,foundnative,repr(replica.chain.contactstate()), repr(replica.chain.vec))

            #print '%-8s %-12s %-12s %-12s'%('replica','foundnative','contact state') 
            #for replica in replicas:
            #    print '%-8d %-12d %s %s'%(replica.repnum,foundnative,repr(replica.chain.contactstate()))

        # Continue to the next production cycle!
        prodstep = prodstep + 1

            
    # write the acceptance ratios to the <rundir>/data directory
    faccept = open(config.DATADIR+'/acceptance','w')
    for i in range(0,len(config.REPLICATEMPS)):
        tmp = str(config.REPLICATEMPS[i]) + '\t'
        tmp = tmp +  str(acceptedsteps[i]) + '\t'
        tmp = tmp +  str(viablesteps[i]) + '\t'
        tmp = tmp +  str(steps[i]) + '\n'


    #Output the status of the simulation one last time...
    print prodstep,'production steps'
    print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s '%('replica','viablesteps','steps','MOVEaccept','viableswaps','swaps','SWAPaccept')
    for rep in range(0,len(steps)):
        print '%-12d %-12d %-12d %-12s %-12d %-12d %-12s '%(rep,viablesteps[rep],steps[rep],'%1.3f'%moveacceptance[rep],swaps[rep],viableswaps[rep],'%1.3f'%swap_acceptance[rep])
    if config.STOPATNATIVE == 1:
        print 'NATIVE CLIST:', nativeclist
    print '%-8s %-12s %-12s %-12s'%('replica','foundnative','contact state', 'chainvec')    
    for replica in replicas:
        print '%-8d %-12d %-12s %12s'%(replica.repnum,foundnative,repr(replica.chain.contactstate()), repr(replica.chain.vec))

    faccept.write(tmp)  
    faccept.close()


    
    # write the last of the trj and ene buffers
    # and close all the open trajectory file handles
    traj.cleanup(replicas)      

