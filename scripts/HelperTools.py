# Functions

import os, sys

import numpy as np
from msmbuilder.MSMLib import *

import random
import string
import math
import os, sys, copy, pickle

import scipy.sparse
import scipy.linalg
import scipy


from scipy.io import mmread, mmwrite

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

def build_MJ_dict():
        """Build a dictionary of Miyazawa contact energies {resi: {resj: energy}} where
        the units are in kT (JMB 1996)."""

        MJtxt = """*** Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro 
    Cys -5.44 -4.99 -5.80 -5.50 -5.83 -4.96 -4.95 -4.16 -3.57 -3.16 -3.11 -2.86 -2.59 -2.85 -2.41 -2.27 -3.60 -2.57 -1.95 -3.07 
    Met 0.46 -5.46 -6.56 -6.02 -6.41 -5.32 -5.55 -4.91 -3.94 -3.39 -3.51 -3.03 -2.95 -3.30 -2.57 -2.89 -3.98 -3.12 -2.48 -3.45 
    Phe 0.54 -0.20 -7.26 -6.84 -7.28 -6.29 -6.16 -5.66 -4.81 -4.13 -4.28 -4.02 -3.75 -4.10 -3.48 -3.56 -4.77 -3.98 -3.36 -4.25 
    Ile 0.49 -0.01 0.06 -6.54 -7.04 -6.05 -5.78 -5.25 -4.58 -3.78 -4.03 -3.52 -3.24 -3.67 -3.17 -3.27 -4.14 -3.63 -3.01 -3.76 
    Leu 0.57 0.01 0.03 -0.08 -7.37 -6.48 -6.14 -5.67 -4.91 -4.16 -4.34 -3.92 -3.74 -4.04 -3.40 -3.59 -4.54 -4.03 -3.37 -4.20 
    Val 0.52 0.18 0.10 -0.01 -0.04 -5.52 -5.18 -4.62 -4.04 -3.38 -3.46 -3.05 -2.83 -3.07 -2.48 -2.67 -3.58 -3.07 -2.49 -3.32 
    Trp 0.30 -0.29 0.00 0.02 0.08 0.11 -5.06 -4.66 -3.82 -3.42 -3.22 -2.99 -3.07 -3.11 -2.84 -2.99 -3.98 -3.41 -2.69 -3.73 
    Tyr 0.64 -0.10 0.05 0.11 0.10 0.23 -0.04 -4.17 -3.36 -3.01 -3.01 -2.78 -2.76 -2.97 -2.76 -2.79 -3.52 -3.16 -2.60 -3.19 
    Ala 0.51 0.15 0.17 0.05 0.13 0.08 0.07 0.09 -2.72 -2.31 -2.32 -2.01 -1.84 -1.89 -1.70 -1.51 -2.41 -1.83 -1.31 -2.03 
    Gly 0.68 0.46 0.62 0.62 0.65 0.51 0.24 0.20 0.18 -2.24 -2.08 -1.82 -1.74 -1.66 -1.59 -1.22 -2.15 -1.72 -1.15 -1.87 
    Thr 0.67 0.28 0.41 0.30 0.40 0.36 0.37 0.13 0.10 0.10 -2.12 -1.96 -1.88 -1.90 -1.80 -1.74 -2.42 -1.90 -1.31 -1.90 
    Ser 0.69 0.53 0.44 0.59 0.60 0.55 0.38 0.14 0.18 0.14 -0.06 -1.67 -1.58 -1.49 -1.63 -1.48 -2.11 -1.62 -1.05 -1.57 
    Asn 0.97 0.62 0.72 0.87 0.79 0.77 0.30 0.17 0.36 0.22 0.02 0.10 -1.68 -1.71 -1.68 -1.51 -2.08 -1.64 -1.21 -1.53 
    Gln 0.64 0.20 0.30 0.37 0.42 0.46 0.19 -0.12 0.24 0.24 -0.08 0.11 -0.10 -1.54 -1.46 -1.42 -1.98 -1.80 -1.29 -1.73 
    Asp 0.91 0.77 0.75 0.71 0.89 0.89 0.30 -0.07 0.26 0.13 -0.14 -0.19 -0.24 -0.09 -1.21 -1.02 -2.32 -2.29 -1.68 -1.33 
    Glu 0.91 0.30 0.52 0.46 0.55 0.55 0.00 -0.25 0.30 0.36 -0.22 -0.19 -0.21 -0.19 0.05 -0.91 -2.15 -2.27 -1.80 -1.26 
    His 0.65 0.28 0.39 0.66 0.67 0.70 0.08 0.09 0.47 0.50 0.16 0.26 0.29 0.31 -0.19 -0.16 -3.05 -2.16 -1.35 -2.25 
    Arg 0.93 0.38 0.42 0.41 0.43 0.47 -0.11 -0.30 0.30 0.18 -0.07 -0.01 -0.02 -0.26 -0.91 -1.04 0.14 -1.55 -0.59 -1.70 
    Lys 0.83 0.31 0.33 0.32 0.37 0.33 -0.10 -0.46 0.11 0.03 -0.19 -0.15 -0.30 -0.46 -1.01 -1.28 0.23 0.24 -0.12 -0.97 
    Pro 0.53 0.16 0.25 0.39 0.35 0.31 -0.33 -0.23 0.20 0.13 0.04 0.14 0.18 -0.08 0.14 0.07 0.15 -0.05 -0.04 -1.75 """

        lines = MJtxt.split('\n')

        mjdict = {}
        residues = lines.pop(0).split()[1:]
        for res in residues:
            mjdict[res] = {}

        for line in lines:
            fields = line.split()
            if len(fields) > 2:
                resi = fields[0]
                for j in range(residues.index(resi),len(residues)):
                    resj = residues[j]
                    mjdict[resi][resj] = float(fields[1+j])
                    mjdict[resj][resi] = mjdict[resi][resj]

        return mjdict


def energy(contact_state, sequence, MJdict):
    """returns the MJ energy of the contacts"""
    ene = 0. 
    for c in contact_state:
        ene += MJdict[sequence[c[0]]][sequence[c[1]]]    
    return ene

def tlc2olc(resname):
    """Return a one-letter code for the amino acid."""

    residues = "Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro".split()
    rescode = "C M F I L V W Y A G T S N Q D E H R K P".split()
    return rescode[ residues.index(resname) ]

def olc2tlc(resname):
    """Return a three-letter code for the amino acid."""

    residues = "Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro".split()
    rescode = "C M F I L V W Y A G T S N Q D E H R K P".split()
    return residues[ rescode.index(resname) ]


def getTimescalesFromTransitionMatrix(T, NumImpliedTimes = 10, LagTime = 1., Verbose=False):
    if Verbose: print('Computing Implied Timescales....')
    EigAns=GetEigenvectors(T,NumImpliedTimes+1,Epsilon=1) #TJL: set Epsilon high, should not raise err here     

    # make sure to leave off equilibrium distribution
    lagTimes = LagTime*np.ones((NumImpliedTimes))
    impTimes = -lagTimes/np.log(EigAns[0][1:NumImpliedTimes+1])

    # save intermediate result in case of failure
    result = np.zeros((NumImpliedTimes, 2))
    result[:,0] = lagTimes
    result[:,1] = impTimes

    return EigAns, result



def draw_index(probs, n_picks=1, UseFastMethod=True):
    """Draw a number (or many numbers, controlled by n_picks), weighted by the probabilities probs."""

    if UseFastMethod:

        t = np.cumsum(probs)
        s = sum(probs)
        return np.searchsorted(t,np.random.rand(n_picks)*s)

    else:
      try:
        nprobs = len(probs)
      except:
        nprobs = probs.shape[0]
      r = np.random.random()
      i = 0
      pthresh = 0.0
      while i < nprobs:
        pthresh += probs[i]
        if r < pthresh:
           return i
        i += 1
      return min(nprobs-1, i)



def reweightTransitionMatrix(T, sequence, MJdict, ContactStateIndices, microContactStates, beta, Verbose=True):
    if Verbose: print('reweighting the transition matrix to reflect the sequence...')
    newT = T.tolil()
    rows, cols = newT.nonzero()
    for k in range(len(rows)):
        if k%10000 == 0:
           if Verbose:  print(k, 'of', len(rows), 'nonzero entries')
        i, j = rows[k], cols[k]
        if ContactStateIndices[i] != ContactStateIndices[j]:
            u_i = energy(microContactStates[i], sequence, MJdict)
            u_j = energy(microContactStates[j], sequence, MJdict)
            newT[i,j] = newT[i,j]*min(1., np.exp(beta*(u_i-u_j)))

    return EstimateTransitionMatrix(newT).tolil()


class MicrostateInfo(object):
    """A class to store microstate info."""
    
    def __init__(self, microstatesFn, NativeContactState=[(0, 11), (2, 11), (4, 9), (4, 11), (6, 9)]):
        """Retrieve and store information about the Microstates from file."""
        
        self.microstates = []
        self.microNumContacts = []
        self.microContactStates = []
        
        self.NativeContactState = NativeContactState

        if os.path.exists(microstatesFn):
            fin = open(microstatesFn, 'r')
            lines = fin.readlines()
            fin.close()
            for line in lines:
                """Example:
        29      [0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 3]       1       [(6, 11)]
        30      [0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0]       1       [(6, 9)]
                """
                fields = line.split('\t') 
                self.microstates.append( eval(fields[1] ))
                self.microNumContacts.append( int(fields[2]))
                self.microContactStates.append( eval(fields[3]))
        
            # compile a list of the 72 unique contact states
            self.uniqueContactStates = []
            for c in self.microContactStates:
                if self.uniqueContactStates.count(c) == 0:
                    self.uniqueContactStates.append(c)
            self.ContactStateIndices = [self.uniqueContactStates.index(self.microContactStates[i]) for i in range(len(self.microstates))]
            self.NativeMicrostateIndex = self.microContactStates.index( self.NativeContactState )
            self.NativeContactStateIndex = self.uniqueContactStates.index(self.NativeContactState)
            print('NativeContactState', self.NativeContactState)
            print('NativeMicrostateIndex', self.NativeMicrostateIndex)
                 
        else:
            print("Can't find file:", microstatesFn, '...exiting')
            raise Exception
            
