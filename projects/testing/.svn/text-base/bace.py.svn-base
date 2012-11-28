#!/usr/bin/env python

import numpy as np
import os, sys

import matplotlib
from matplotlib import pyplot as plt

def draw_index(probs):
    """Draw a number, weighted by the probabilities probs."""

    r = np.random.random()
    i = 0
    pthresh = 0.0
    while i < len(probs):
        pthresh += probs[i]
        if r < pthresh:
           return i
        i += 1
    return min(len(probs)-1, i)

def logPdiff(C1, C2, verbose=False):
    """Return the BACE ln(Pdiff/Psame) or two count matrices."""

    Csum = C1+C2
    result = 0.
    for i in range(C1.shape[0]):
        Ci, Cj, Cij = C1[i,:].sum(), C2[i,:].sum(), Csum[i,:].sum()
        pi, pj, q = C1[i,:]/Ci, C2[i,:]/Cj, Csum[i,:]/Cij
        if i == 0 and verbose:
            print '\tpi', pi,
            print '\tpj', pj,
            print '\tq', q,
            print '\trelativeEntropy(pi,q)', relativeEntropy(pi,q)
            print '\trelativeEntropy(pj,q)', relativeEntropy(pj,q)
        result += Ci*relativeEntropy(pi,q) + Cj*relativeEntropy(pj,q)
    return result
   
def relativeEntropy(p,q):
    """Returns the relative entropy of two distributions, p and q.
    NOTE: This function assumes p and q are properly normalized."""

    Dterms = p*np.log(p) - p*np.log(q)
    return Dterms[~np.isnan(Dterms)].sum()


 

N = 10
T1 = np.random.random( (N,N) ) 
T2 = T1
#T2 = 0.999*T1 + 0.001*np.random.random( (N,N) ) 
for i in range(N):
    T1[i,:] = T1[i,:]/T1[i,:].sum()
    T2[i,:] = T2[i,:]/T2[i,:].sum()
print T1, T2

C1 = np.zeros( T1.shape )
C2 = np.zeros( T2.shape )

# draw m counts as a Markov Chain
total = 0
m = 1000
steps = range(100000)
data = []
for step in steps: 

  i = np.random.randint(N)
  for k in range(m):
      j = draw_index(T1[i,:])
      C1[i,j] += 1.
      i = j

  i = np.random.randint(N)
  for k in range(m):
      j = draw_index(T2[i,:])
      C2[i,j] += 1.
      i = j

  total += m
  #print 'C1', C1
  #print 'C2', C2
  logP = logPdiff(C1, C2)
  print total, 'logPdiff(C1, C2)', logP
  data.append(logP)

plt.figure()
plt.plot(steps, data)
plt.show()


