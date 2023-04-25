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


# Functions

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


def getTimescales(Counts, NumImpliedTimes = 20, LagTime = 1., Verbose=False):
    if Verbose: print 'Computing Implied Timescales....'
    newT = EstimateTransitionMatrix(Counts)
    EigAns=GetEigenvectors(newT,NumImpliedTimes+1,Epsilon=1) #TJL: set Epsilon high, should not raise err here     

    # make sure to leave off equilibrium distribution
    lagTimes = LagTime*np.ones((NumImpliedTimes))
    impTimes = -lagTimes/np.log(EigAns[0][1:NumImpliedTimes+1])

    # save intermediate result in case of failure
    result = np.zeros((NumImpliedTimes, 2))
    result[:,0] = lagTimes
    result[:,1] = impTimes

    return EigAns, result


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

def reweightTransitionMatrix(T, sequence, MJdict, ContactStateIndices, microContactStates, beta, Verbose=True):
    if Verbose: print 'reweighting the transition matrix to reflect the sequence...'
    newT = T.tolil()
    rows, cols = newT.nonzero()
    for k in range(len(rows)):
        if k%10000 == 0:
           if Verbose:  print k, 'of', len(rows), 'nonzero entries'
        i, j = rows[k], cols[k]
        if ContactStateIndices[i] != ContactStateIndices[j]:
            u_i = energy(microContactStates[i], sequence, MJdict)
            u_j = energy(microContactStates[j], sequence, MJdict)
            newT[i,j] = newT[i,j]*min(1., np.exp(beta*(u_i-u_j)))

    return EstimateTransitionMatrix(newT).tolil()



################
# Main program

# Read in transition matix
microTFn = 'microT.mtx'
if os.path.exists(microTFn):
    T = mmread(microTFn)
    print T
else:
    print "Can't find file:", microTFn, '...exiting'
    sys.exit(1)

# Read in microstate info
microstatesFn = 'microstates.dat'
microstates = []         # list of vecs
microContactStates = [] # list of contact states for each 
microNumContacts  = []   # list of number of contacts

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
        microstates.append( eval(fields[1] ))
        microNumContacts.append( int(fields[2]))
        microContactStates.append( eval(fields[3]))

    # compile a list of the 72 unique contact states
    uniqueContactStates = []
    for c in microContactStates:
        if uniqueContactStates.count(c) == 0:
            uniqueContactStates.append(c)
    ContactStateIndices = [uniqueContactStates.index(microContactStates[i]) for i in range(len(microstates))]
    NativeContactState = [(0, 11), (2, 11), (4, 9), (4, 11), (6, 9)]
    NativeMicrostateIndex = microContactStates.index( NativeContactState )
    NativeContactStateIndex = uniqueContactStates.index(NativeContactState)
    print 'NativeContactState', NativeContactState
    print 'NativeMicrostateIndex', NativeMicrostateIndex
         
else:
    print "Can't find file:", microstatesFn, '...exiting'
    sys.exit(1)


############################

LagTime = 1

if (1):
    print '#stability\tslowest\tnextslowest\tslowest2\tnextslowest2'
    ########

    if os.path.exists('newT1.mtx'):
        print 'Reading from newT1.mtx ...'
        newT1 = mmread('newT1.mtx')
    else:

        # Build sequence information
        MJdict = build_MJ_dict()
        residues = "Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro".split()
        resIndices = [0,2,4,6,9,11]
        sequence = {}
        initial_sequence = ['Ala', 'Ala', 'Ala', 'Ala', 'Tyr', 'Ala']
        for i in resIndices:
            sequence[i] = initial_sequence[resIndices.index(i)]
        print 'sequence = ', sequence
        beta = 1.0
        newT1 = reweightTransitionMatrix(T, sequence, MJdict, ContactStateIndices, microContactStates, beta)

        for i in range(LagTime-1):
            newT1 = newT1.dot( newT1 ) 

        mmwrite('newT1.mtx', newT1, comment='')

    try: 
        # get the Implied Timescales
        print 'Calculating timescales for the unperturbed T...'
        EigAns, result = getTimescales(newT1, LagTime=LagTime)
        EnoughCounts = True
    except:
        EnoughCounts = False
     
    if EnoughCounts: 
        # calculate the stability of the native state
        INative = np.array(ContactStateIndices)==NativeContactStateIndex
        stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )

        # print stability and timescale info 
        print '%16.8f\t%16.8f\t%16.8f'%(stability, result[0,1], result[1,1])

        lambda2_unperturbed = np.real(EigAns[0][1])
    else:
        print 'not enough counts yet...'

if (1):
    ########
    print 'Perturbing the matrix....'

    if os.path.exists('newT2.mtx'):
        print 'Reading from newT2.mtx ...'
        newT2 = mmread('newT2.mtx')
    else:

        # Build sequence information
        MJdict = build_MJ_dict()
        residues = "Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro".split()
        resIndices = [0,2,4,6,9,11]
        sequence = {}
        initial_sequence = ['Ala', 'Ala', 'Ala', 'Ala', 'Ala', 'Ala']
        for i in resIndices:
            sequence[i] = initial_sequence[resIndices.index(i)]
        print 'sequence = ', sequence

        # Reweight the count matrix to reflect our chosen sequence
        # Recall the microT.mtx (T) is just ones and zeros (i.e. counts)
        beta = 1.0
        newT2 = reweightTransitionMatrix(T, sequence, MJdict, ContactStateIndices, microContactStates, beta) 

        for i in range(LagTime-1):
            newT2 = newT1.dot( newT2 ) 

        mmwrite('newT2.mtx', newT2, comment='')

    try:
        # get the Implied Timescales
        print 'Calculating timescales for the perturbed T ...'
        newEigAns, newresult = getTimescales(newT2, LagTime=LagTime)
        EnoughCounts = True
    except:
        EnoughCounts = False

    if EnoughCounts:
        # calculate the stability of the native state
        INative = np.array(ContactStateIndices)==NativeContactStateIndex
        stability = np.real( newEigAns[1][INative,0]/newEigAns[1][:,0].sum() )

        # print stability and timescale info 
        print '%16.8f\t%16.8f\t%16.8f'%(stability, newresult[0,1], newresult[1,1])

        lambda2_perturbed = np.real(newEigAns[0][1])
    else:
        print 'not enough counts yet...' 


if (1):
    # NOW, let's try perturbation theory
    UseSymmetrizedT = True 

    if UseSymmetrizedT:
        e1 = EigAns[1][:,0]
        equilpops = e1/e1.sum()
        #Peq1 = scipy.sparse.dia_matrix((equilpops**0.5,0),newT1.shape).tocsr()
        Peq1 = scipy.sparse.spdiags( [equilpops**0.5], [0], newT1.shape[0], newT1.shape[1], format="csr")
        print 'Peq1', Peq1
        #Peq2 = scipy.sparse.dia_matrix((equilpops**-0.5,0),newT1.shape).tocsr()
        Peq2 = scipy.sparse.spdiags( [equilpops**-0.5], [0], newT1.shape[0], newT1.shape[1], format="csr")
        print 'Peq2', Peq2
        Tsym = Peq1.dot( newT1.dot(Peq2) ).tocsr()
        print 'Tsym', Tsym

        #Tsym = np.real((Tsym + Tsym.transpose())*0.5)  # naive symmetrization (because small numbers are off)
        #print 'Tsym', Tsym

        n=Tsym.shape[0]
        NumImpliedTimes = 20
        NumEig = NumImpliedTimes + 1
        if NumEig>n:
            raise Exception("You cannot calculate %d Eigenvectors from a %d x %d matrix"%(NumEig,n,n))

        from scipy.sparse.linalg import eigsh
        print 'computing eigSolution...'
        if scipy.sparse.issparse(Tsym):
            eigSolution = eigsh(Tsym, k=10, which='LR')
            #eigSolution = sparseEigenSym(Tsym.tocsr(), NumEig, which="LR", maxiter=1000000, tol=0)
            #eigSolution = sparseEigen(Tsym.tocsr(), NumEig, which="LR", maxiter=1000000, tol=0)
        else:
            eigSolution=eig(Tsym)
        print 'Ordering...'
        Ord=np.argsort(-np.real(eigSolution[0]))

        elambda = eigSolution[0][Ord]
        eV=eigSolution[1][:,Ord]

        print 'check that the eV are properly normalized...'
        print 'np.dot(eV[:,0],eV[:,0])', np.dot(eV[:,0],eV[:,0])
        print 'np.dot(eV[:,1],eV[:,0])', np.dot(eV[:,1],eV[:,0])
        print 'np.dot(eV[:,1],eV[:,1])', np.dot(eV[:,1],eV[:,1])


        dT = newT2 - newT1
        e2 = eV[:,1]
        #dlambda = LeftE2.transpose().dot(  dT.dot(e2) )
        dlambda = e2.transpose().dot(  dT.dot(e2) )


     
    else:
        n=newT1.shape[0]
        NumImpliedTimes = 20
        NumEig = NumImpliedTimes + 1
        if NumEig>n:
            raise Exception("You cannot calculate %d Eigenvectors from a %d x %d matrix"%(NumEig,n,n))
        if scipy.sparse.issparse(newT1):
            eigSolution = sparseEigen(newT1.tocsr(), NumEig, which="LR", maxiter=1000000, tol=0)
            LeftEigSolution = sparseEigen(newT1.transpose().tocsr(), NumEig, which="LR", maxiter=1000000)
        else:
            eigSolution=eig(T)
            LeftEigSolution = eig(T.transpose())
        Ord=np.argsort(-np.real(eigSolution[0]))
        LeftOrd=np.argsort(-np.real(LeftEigSolution[0]))


        print 'Check to see that left and right evalues are the same:'
        elambda = eigSolution[0][Ord]
        LeftElambda = LeftEigSolution[0][LeftOrd]
        print 'elambda', elambda, 'LeftElambda', LeftElambda

        eV=eigSolution[1][:,Ord]
        LeftEV = LeftEigSolution[1][:,LeftOrd]

        print 'check that the eV are properly normalized...'
        print 'np.dot(eV[:,1],eV[:,1])', np.dot(eV[:,1],eV[:,1])
        print 'np.dot(LeftEV[:,1],LeftEV[:,1])', np.dot(LeftEV[:,1],LeftEV[:,1])
        print 'np.dot(LeftEV[:,1],eV[:,1])', np.dot(LeftEV[:,1],eV[:,1])


        dT = newT2 - newT1
        e2 = eV[:,1]
        LeftE2 = LeftEV[:,1]
        dlambda = np.real(LeftE2.transpose().dot(  dT.dot(e2) ))

    print 'lambda2_unperturbed', lambda2_unperturbed, 'lambda2_perturbed', lambda2_perturbed
    print 'diff = ', lambda2_perturbed-lambda2_unperturbed
    print 'Predicted by PT:', dlambda
    print '--------------'
    LagTime = 1.0
    tau2_unperturbed = -LagTime/np.log(lambda2_unperturbed)
    tau2_perturbed = -LagTime/np.log(lambda2_perturbed)
    print 'tau_2 unperturbed =', tau2_unperturbed, 'tau_2 perturbed =', tau2_perturbed
    print 'diff = ', tau2_perturbed - tau2_unperturbed  
    print 'Predicted by PT:', -LagTime/np.log(lambda2_unperturbed + dlambda) + LagTime/np.log(lambda2_unperturbed)

       
