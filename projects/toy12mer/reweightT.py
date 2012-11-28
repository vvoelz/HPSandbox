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


usage = """Usage: reweightT.py infile.mtx outfile.mtx microstates.dat sequence 

    sequence - a one-letter amino acid sequence, for the six hydrophobic(H) residues
               in the HPHPHPHPPHPH sequence:  AAAAAA  for example

    Will write files [sequence].dat containing the stability and slowest, next-slowest implied timescales
    """


sys.path.append('../../scripts')
from HelperTools  import *

################
# Main program

if len(sys.argv) < 5:
    print usage
    sys.exit(1)

microTFn            = sys.argv[1]
reweighted_microTFn = sys.argv[2]
microstatesFn       = sys.argv[3]
sequence            = sys.argv[4]
outfile = sequence + '.dat'

# Read in transition matrix
if os.path.exists(microTFn):
    T = mmread(microTFn)
    print T
else:
    print "Can't find file:", microTFn, '...exiting'
    sys.exit(1)

# Read in microstate info
microstates = []         # list of vecs
microContactStates = [] # list of contact states for each 
microNumContacts  = []   # list of number of contacts

m =  MicrostateInfo(microstatesFn)

# Build sequence information
if os.path.exists(reweighted_microTFn):
        newT = mmread(reweighted_microTFn)
        print newT
else:
        MJdict = build_MJ_dict()
        residues = "Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro".split()
        resIndices = [0,2,4,6,9,11]
        seqdict = {} # {resindex: 'Cys'} e.g.
        initial_sequence = [olc2tlc(sequence[i]) for i in range(len(sequence))]
        for i in resIndices:
            seqdict[i] = initial_sequence[resIndices.index(i)]
        print 'seqdict = ', seqdict
        beta = 1.0
        newT = reweightTransitionMatrix(T, seqdict, MJdict, m.ContactStateIndices, m.microContactStates, beta)

        mmwrite(reweighted_microTFn, newT, comment=sequence)

CalculateOriginalTimescales = True #False

if CalculateOriginalTimescales:
    print 'Original T:'
    EigAns, result = getTimescalesFromTransitionMatrix(T, NumImpliedTimes = 10)
    print EigAns
    print result

    # calculate the stability of the native state
    INative = np.array(m.ContactStateIndices)==m.NativeContactStateIndex
    stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )
    # print stability and timescale info 
    print '#stability\tslowest\tnext-slowest'
    print '%8.3f\t%16.8f\t%16.8f'%(stability, result[0,1], result[1,1])


print 'reweighted T:'
EigAns, result = getTimescalesFromTransitionMatrix(newT, NumImpliedTimes = 10)
print EigAns
print result

# calculate the stability of the native state
INative = np.array(m.ContactStateIndices)==m.NativeContactStateIndex
stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )
# print stability and timescale info 
header = '#sequence\tstability\tslowest\tnext-slowest'
line = '%s\t%8.3f\t%16.8f\t%16.8f'%(sequence, stability, result[0,1], result[1,1])
print header+'\n'+line
# write stability and timescale info to output file
fout = open(outfile,'w')
fout.write(header+'\n')
fout.write(line+'\n')
fout.close()

# Calculate the timescale and stability of the macrostate transition matrix

NumContactStates = len(m.uniqueContactStates)
C_macro  = scipy.sparse.lil_matrix((int(NumContactStates),int(NumContactStates)))

newT_micro = newT.tolil()
print 'newT_micro.rows', newT_micro.rows
print 'newT_micro.data', newT_micro.data

for i in range(newT.shape[0]):
    for j in range(len(newT_micro.rows[i])):
        k = newT_micro.rows[i][j]
        C_macro[m.ContactStateIndices[i], m.ContactStateIndices[k]] += newT_micro.data[i][j]

T_macro = EstimateTransitionMatrix(C_macro)
print 'T_macro', T_macro

mmwrite(reweighted_microTFn.replace('micro','macro'), T_macro, comment='macrostate '+sequence)

print 'reweighted T_macro:'
EigAns, result = getTimescalesFromTransitionMatrix(T_macro, NumImpliedTimes = 69)
print EigAns
print result

# calculate the stability of the native state from the MACROSTATE 
INative = np.array(range(int(NumContactStates)))==m.NativeContactStateIndex
stability = np.real( EigAns[1][INative,0]/EigAns[1][:,0].sum() )
# print stability and timescale info 
print '#stability\tslowest\tnext-slowest'
print '%8.3f\t%16.8f\t%16.8f'%(stability, result[0,1], result[1,1])

