#! /usr/bin/python
#
# Monty.py      
#    
# PURPOSE:
#       An object to implement Monte Carlo movesets on a lattice chain 
#       
# NOTES:
#       Fixed random number generator problem! don't use g.random(), as seed will be the same...        

from Config import *
from Chain import *
from Replica import *
from Trajectory import *

import random
import string
import math
import sys
import os

randseed = 422

#########################################################################

class Monty:
    """A collection of functions to perform Monte Carlo move-set operations on an HP lattice Chain object."""

    def __init__(self,config,temp,chain):
        """Initialize the Monte Carlo object..."""
    
	print '\tcreating Monty.py object....'
	self.g = random.Random(config.randseed)               # NOTE:  random.Random is an object that
                                                              # needs an integer seed for initialization
    
	self.movesets = ['MC1','MC2','MC3','MC4']             # The names of the available Monte Carlo movesets
	self.temp = temp                                      # The temperature (in K)
	self.tempfromrep = config.REPLICATEMPS.index(temp)    # The replica number with this temperature
	self.epsilon = config.epsilon                         # The energetic strength of a contact 
	self.k = config.k				      # Boltzmann's constant
							      # (both of these are copies from Config() )
	self.restraint = DistRestraint(config.RESTRAINED_STATE,config.KSPRING)
	self.lastenergy = self.energy(chain) + self.restraint.energy(chain)
	
	
    def move1(self,replica):
	"""Apply moveset 'MC1' to the chain:
	(i)  three-bead flips
	(ii) end flips
	
	REFERENCE: Dill and Chan, 1994, 1996.
"""
    
	r = random.random()
	s = random.random()
     
	vecindex = int(math.floor((replica.chain.n - 1.0001)*r))
	
	if s < 0.5:
	    direction = 1      
	else:
	    direction = -1      
 
	### if the vec index is on the end, do an end flip
	if (vecindex == 0) | (vecindex == len(replica.chain.nextvec)-1):
	    replica.chain.nextvec[vecindex] = (replica.chain.nextvec[vecindex] + direction) % 4

	### else, do a three-bead flip, i.e. switch two adjacent {n,e,w,s} directions
	else:
	    tmp = replica.chain.nextvec[vecindex] 
	    replica.chain.nextvec[vecindex] = replica.chain.nextvec[vecindex+1] 
	    replica.chain.nextvec[vecindex + 1] = tmp
	
	replica.chain.nextcoords = replica.chain.vec2coords(replica.chain.nextvec)
	replica.chain.nextviable = replica.chain.viability(replica.chain.nextcoords)


    def move2(self,replica):
	"""Apply moveset MC2 to the chain:
	(i)   three-bead flips
	(ii)  end flips
	(iii) crankshaft moves
	(iv)  rigid rotations
	
	REFERENCE:  Dill and Chan, 1994, 1996
"""
    
	r = random.random()
	s = random.random()
	t = random.random()
	     
	vecindex = int(math.floor((replica.chain.n - 1.0001)*r))
	
	if s < 0.5:
	    direction = 1      
	else:
	    direction = -1      

        # if possible, 1/3 of the time do a three-bead flip (dirs must be different)
	if (t < 0.33333) & (vecindex < len(replica.chain.nextvec)-1):
	  tmp1 = replica.chain.nextvec[vecindex] 
	  tmp2 = replica.chain.nextvec[vecindex+1]
	  if (tmp1 != tmp2):
	      replica.chain.nextvec[vecindex] = tmp2 
	      replica.chain.nextvec[vecindex + 1] = tmp1
	  else: 
	      ### default: do a rigid rotation
	      for v in range(vecindex,len(replica.chain.nextvec)):
		replica.chain.nextvec[v] = (replica.chain.nextvec[v] + direction) % 4
		# print 'chain.nextvec',chain.nextvec

	# if possible, 1/3 of the time do a crankshft (1st and 3rd dirs must be different)
	elif (t < 0.66666) & (vecindex < len(replica.chain.nextvec)-2):
	  tmp1 = replica.chain.nextvec[vecindex] 
	  tmp2 = replica.chain.nextvec[vecindex+2]
	  if (t < 0.66666) & (tmp1 != tmp2):
	      ### crankshaft move
	      replica.chain.nextvec[vecindex] = tmp2 
	      replica.chain.nextvec[vecindex + 2] = tmp1
	  else: 
	      ### default: do a rigid rotation
	      for v in range(vecindex,len(replica.chain.nextvec)):
		replica.chain.nextvec[v] = (replica.chain.nextvec[v] + direction) % 4
		# print 'chain.nextvec',chain.nextvec

	else: 
	    ### default: do a rigid rotation
	    for v in range(vecindex,len(replica.chain.nextvec)):
	      replica.chain.nextvec[v] = (replica.chain.nextvec[v] + direction) % 4
	      # print 'chain.nextvec',chain.nextvec
	
	replica.chain.nextcoords = replica.chain.vec2coords(replica.chain.nextvec)
	replica.chain.nextviable = replica.chain.viability(replica.chain.nextcoords)


    def move3(self,replica):
	"""Apply moveset 'MC3' to the chain.
        This is just a simple set to change the direction of a single chain link.
        Example:
            [0,0,0,0,0] --> [0,0,1,0,0]
        where {0,1,2,3}={n,e,s,w} direction

	About 5% viable moves are expected."""
    
	r = self.random.random()
	s = self.random.random()
     
	vecindex = int(math.floor((replica.chain.n - 1.0001)*r))
	
	if s < 0.5:
	    direction = 1      
	else:
	    direction = -1      
  
	replica.chain.nextvec[vecindex] = (replica.chain.nextvec[vecindex] + direction) % 4
	replica.chain.nextcoords = replica.chain.vec2coords(replica.chain.nextvec)
	replica.chain.nextviable = replica.chain.viability(replica.chain.nextcoords)


    def move4(self,replica):
	"""Apply moveset 'MC4' to the chain:
	This is another vert simple moveset, to just change one angle in a rigid rotation
	Like 'MS3', this generates about 5% viable moves."""
    
	r = self.g.random()
	s = self.g.random()
     
	vecindex = int(math.floor((replica.chain.n - 1.0001)*r))
	
	if s < 0.5:
	    direction = 1      
	else:
	    direction = -1      
  
	### a rigid rotation
	for v in range(vecindex,len(replica.chain.nextvec)):
	    replica.chain.nextvec[v] = (replica.chain.nextvec[v] + direction) % 4
	
	replica.chain.nextcoords = replica.chain.vec2coords(replica.chain.nextvec)
	replica.chain.nextviable = replica.chain.viability(replica.chain.nextcoords)



    def metropolis(self,replica):
	"""Accept Chain.nextvec over Chain.vec according to a Metropolis criterion."""
    
	randnum = self.g.random()
	
	# accept with Metroplis criterion      
	thisenergy = self.energy(replica.chain) + self.restraint.energy(replica.chain)  
	boltzfactor = math.exp((thisenergy - self.lastenergy)/(self.k*self.temp))
	
	if randnum < boltzfactor:
	    # update the chain
	    for i in range(0,len(replica.chain.vec)):
		replica.chain.vec[i] = replica.chain.nextvec[i]
	    for i in range(0,len(replica.chain.coords)):
		replica.chain.coords[i] = replica.chain.nextcoords[i]
            replica.chain.viable = replica.chain.nextviable
	    # update the lastenergy
	    self.lastenergy = thisenergy
	    return 1
	else:
	    return 0
     
    def energy(self,chain):
	"""Calculate potential energy of the chain."""
	
	num = 0.0
	for c in range(0,len(chain.coords)-1):
	    for d in range((c+3),len(chain.coords)):
	      if chain.hpstring[c] == 'H':
		if chain.hpstring[d] == 'H':
		  if (abs(chain.coords[c][0]-chain.coords[d][0]) + abs(chain.coords[c][1]-chain.coords[d][1])) == 1:
		    num = num + 1.0
	return num*self.epsilon



#########################################################################


class DistRestraint:
    """ For now, this is a harmonic constraint over a squared distance D = d^2
     where D = sum_{i,j} d^2_ij over all contacts."""
    
    def __init__(self,contacts,kspring):
        """Initialize the DistRestraint object"""
	self.contacts = contacts        # a list of duples (start of chain is 0)
	self.kspring = kspring          # spring constant for restraint
					# (J/[lattice space]^2)

    def energy(self,chain):
        """ return the energy of the distance restraint"""
	return self.kspring*self.D(chain)
	
    def D(self,chain):
	"""Return the sum of squared-distances over the selected contacts."""
	D = 0.0
	for i in range(0,len(self.contacts)):
	    c = self.contacts[i][0]
	    d = self.contacts[i][1]
	    D = D + (chain.coords[c][0]-chain.coords[d][0])*(chain.coords[c][0]-chain.coords[d][0])
	    D = D + (chain.coords[c][1]-chain.coords[d][1])*(chain.coords[c][1]-chain.coords[d][1])

	return D

	
