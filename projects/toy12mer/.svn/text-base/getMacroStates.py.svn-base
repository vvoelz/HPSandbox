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


usage = """Usage: getMacrostates.py microstates.dat 

    """

sys.path.append('../../scripts')
from HelperTools  import *

################
# Main program

if len(sys.argv) < 2:
    print usage
    sys.exit(1)

microstatesFn       = sys.argv[1]

m = MicrostateInfo(microstatesFn)
print m.uniqueContactStates


