import numpy as np

# import dense eignesolvers
#import scipy.linalg
#eig=scipy.linalg.eig
from numpy import linalg as LA
eig = LA.eig

from scipy.io import mmread, mmwrite

# Load sparse eigensolvers
import scipy.sparse.linalg.eigen.arpack as arpack
sparseEigen = arpack.eigs
#sparseEigen = arpack.eigsh  # for symmetric matrices

N = 5 
UseSparse = False
if UseSparse:
    C = scipy.sparse.lil_matrix((N,N),dtype=np.float64)
else:
    C = np.zeros((N,N),dtype=np.float64)

# Fill the sparse matrix with a diagonal, and some random entries
NumEntries = 10*N
for i in range(N):
    C[i,i] = 1.0
for k in range(NumEntries):
    i,j = np.random.randint(N), np.random.randint(N) 
    C[i,j] += np.random.random() 
    C[j,i] += C[i,j]
print 'C', C

# make a row-normalized Transition Matrix
from msmbuilder.MSMLib import *
T = EstimateTransitionMatrix(C)

# symmetrize this matrix 
if UseSparse:
    eigSolution = sparseEigen(T.transpose(), NumEig, which="LM")
    print 'eigSolution', eigSolution
    e1 = eigSolution[1][:,0]
    equilpops = e1/e1.sum()
    Peq1 = scipy.sparse.spdiags( [equilpops**0.5], [0], T.shape[0], T.shape[1], format="csr")
    print 'Peq1', Peq1
    Peq2 = scipy.sparse.spdiags( [equilpops**-0.5], [0], T.shape[0], T.shape[1], format="csr")
    print 'Peq2', Peq2
    Tsym = Peq1.dot( T.dot(Peq2) ).tocsr()
    print 'Tsym', Tsym
else:
    eigSolution = eig(T.transpose(), right=True)
    print 'eigSolution', eigSolution
    e1 = eigSolution[1][:,0]
    equilpops = e1/e1.sum()
    Peq1 = np.diag(equilpops**0.5)
    print 'Peq1', Peq1
    Peq2 = np.diag(equilpops**-0.5)
    print 'Peq2', Peq2
    Tsym = Peq1.dot( T.dot(Peq2) )
    print 'Tsym', Tsym


    print 'equilpops', equilpops
    print '(T.transpose()).dot(equilpops)', (T.transpose()).dot(equilpops)


sys.exit(1)

#from numpy import linalg as LA
print 'LA.cond(T)', LA.cond(T)

print 'T', T
NumEig = min(6,N-2)
UseSparse = False

if UseSparse:
    eigSolution     = sparseEigen(T,             NumEig, which="LM", maxiter=10000000, sigma=1.0, OPpart='r')
    LeftEigSolution = sparseEigen(T.transpose(), NumEig, which="LM", maxiter=10000000, sigma=1.0, OPpart='r')
else:
    # help: eig(a, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False)
    eigSolution     = eig(T, right=True)
    LeftEigSolution = eig(T, left=True)

Ord=np.argsort(-np.real(eigSolution[0]))
LeftOrd=np.argsort(-np.real(LeftEigSolution[0]))

eV=eigSolution[1][:,Ord]
LeftEV = LeftEigSolution[1][:,LeftOrd]



print 'Check to see that left and right evalues are the same:'
elambda = eigSolution[0][Ord]
LeftElambda = LeftEigSolution[0][LeftOrd]
print 'index\telambda\tLeftElambda'
for i in range(NumEig):
    print i, np.real(elambda[i]), np.real(LeftElambda[i])

print 'Checking evectors:'
print 'index\telambda\tLeftElambda'
for i in range(NumEig):
    print i, LeftEV[:,i], eV[:,i]


for k in range(NumEig):
    eV[:,k]  = eV[:,k]/np.dot(eV[:,k],eV[:,k])**0.5
    LeftEV[:,k]  = LeftEV[:,k]/np.dot(LeftEV[:,k],LeftEV[:,k])**0.5
print 'check that the eV are properly normalized...'
print 'np.dot(eV[:,1],eV[:,1])', np.dot(eV[:,1],eV[:,1])
print 'np.dot(LeftEV[:,1],LeftEV[:,1])', np.dot(LeftEV[:,1],LeftEV[:,1])
print 'np.dot(LeftEV[:,1],eV[:,1])', LeftEV[:,1].dot(eV[:,1])
print 'np.dot(LeftEV[:,0],eV[:,1])', np.dot(LeftEV[:,0],eV[:,1])
print 'np.dot(LeftEV[:,1],eV[:,0])', np.dot(LeftEV[:,1],eV[:,0])


