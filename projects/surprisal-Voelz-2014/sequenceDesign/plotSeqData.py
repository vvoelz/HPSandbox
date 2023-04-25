#!/usr/bin/env python

import os, sys
import numpy as np
from matplotlib import pyplot

fin = open('seqdata.dat')
example ="LYYYYY	   0.685	 357287.43417965	  61939.82758913	True"
lines = fin.readlines()
fin.close()

seq, stability, lograte, lograte_gap = [],[],[],[]
for line in lines:
    fields = line.split()
    if len(fields) > 0:
        seq.append(fields[0])
        stability.append( float(fields[1]))
        lograte.append(np.log(float(fields[2])))
        lograte_gap.append(lograte[-1] - np.log(float(fields[3])))

print seq[np.argmax(stability)]


pyplot.figure()
pyplot.subplot(2,2,1)
pyplot.plot(stability, lograte,'k.')
pyplot.subplot(2,2,2)
pyplot.plot(stability, lograte_gap,'k.')
pyplot.show()

