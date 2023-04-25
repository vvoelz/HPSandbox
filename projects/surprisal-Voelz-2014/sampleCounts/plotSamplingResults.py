import os, sys
import numpy as np

from scipy import loadtxt

import matplotlib
from pylab import *


results = loadtxt('uniform.dat')
maxsteps = max(results[:,0])


figure()

subplot(2,2,1)
plot( results[:,0], results[:,1] )
xlabel('number of transition counts')
ylabel('native stability')

subplot(2,2,2)
semilogy( results[:,0], results[:,2:] )
xlabel('number of transition counts')
ylabel('two slowest implied timescales')
axis([0,maxsteps,1e3,1e6])

subplot(2,2,3)
semilogy( results[:,0], results[:,2] )
xlabel('number of transition counts')
ylabel('folding time $\\tau_1$')
axis([0,maxsteps,1e3,1e6])

subplot(2,2,4)
noAdaptive = loadtxt('noAdapt.dat')
yesAdaptive = loadtxt('yesAdapt.dat')
semilogy( noAdaptive[:,0], noAdaptive[:,2] )
hold(True)
semilogy( yesAdaptive[:,0], yesAdaptive[:,2] )
xlabel('number of transition counts')
ylabel('folding time $\\tau_1$')
axis([0,max(noAdaptive[:,0]),0.1,10.])
legend(['noAdaptive', 'yesAdaptive'])

show()






