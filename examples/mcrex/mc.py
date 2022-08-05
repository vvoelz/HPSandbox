#! /usr/bin/python

# mc.py   - a really crappy test program!

from random import Random

import math
import sys
import os

########### Functions ##############

def numcontacts(coords):

    num = 0
    for c in range(0,len(coords)-1):
      for d in range((c+3),len(coords)):

        if (  abs(coords[c][0]-coords[d][0]) + abs(coords[c][1]-coords[d][1])   ) == 1:
	    num = num + 1

    return num
    

########### Program ##################3

g = Random(420)
g.random()

T = 1.0  		# unitless (for now)
k = 1.0	    		# unitless (for now)
epsilon = -1.0*k*T	# units of kT (for now)

Z = []   	# The partition function for each num of contacts
F = []   	# The free energy for each num of contacts
maxcontacts = 10
for i in range(0,maxcontacts+1):
    Z.append(0)    
    F.append(0)    

vec = [0,0,0,0,0,0,0,0,0]	# 10-mer vector angles {0,1,2,3} =  {n,w,s,e}

mcsteps = 10000

for i in range(0,mcsteps):

     r = g.random()
     s = g.random()
     
     vecindex = int(math.floor(8.999*r))
     if s < 0.5:
        direction = 1      
     else:
        direction = -1      

     oldvec = []
     for v in range(0,len(vec)):
         oldvec.append(vec[v])
	 
     vec[vecindex] = ( vec[vecindex] + direction) % 4
    
     coords = [(0,0)]

     x = 0
     y = 0
     for i in range(0,len(vec)):
	 if vec[i] == 0:
	     y = y + 1
	 if vec[i] == 1:
	     x = x + 1
	 if vec[i] == 2:
	     y = y - 1
	 if vec[i] == 3:
	     x = x - 1
  	 coords.append((x,y))


     ### ...and determine whether it's viable	
     viable = 1
     for c in coords:
	 if coords.count(c) > 1:
	     viable = 0
	     break

     if viable == 0:
         for v in range(0,len(vec)):
	     vec[v] = oldvec[v]
	     
     else:
         ### accept with Metroplis criterion	 
	 n = numcontacts(coords)
	 boltzfactor = math.exp(n*epsilon/(k*T))
	 if g.random() > boltzfactor:
	     vec[v] = oldvec[v]
	 else:
	     Z[n] = Z[n] + boltzfactor        
     
     # print 'Z[c] =',Z
     tmp = ''
     for c in range(0,len(Z)):
         F[c] = -1.0*k*T*math.log(max(Z[c],0.0001)) - -1.0*k*T*math.log(max(Z[0],0.0001))
         tmp = tmp + str(F[c]) + '\t' 
     	 
     print tmp
     
