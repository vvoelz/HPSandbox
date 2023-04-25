import os, sys

from matplotlib import pyplot
from scipy import loadtxt
import numpy as np

usage = """Usage: python plot_observable.py obsfile1 obsfile2 .... obsfileN"""

if len(sys.argv) < 2:
    print usage
    sys.exit(1)

obsfiles = sys.argv[1:]
nfiles = len(obsfiles)

pyplot.figure()

for i in range(len(obsfiles)):
    obsfile = obsfiles[i]
    fin = open(obsfile,'r')
    lines = fin.readlines()
    lines.pop(0)
    data = np.array([ [int(field) for field in line.split()[0:3] ] for line in lines ])

    pyplot.subplot(nfiles,1,i+1)
    pyplot.plot(data[:,0], data[:,1])

pyplot.show()


