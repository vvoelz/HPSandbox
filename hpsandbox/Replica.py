#
# Replica.py

#from Config import *
#import Chain
#import Monty
#from Trajectory import *

import random
import string
import math
import sys
import os
from .Chain import Chain
from .Monty import Monty

class Replica:
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self,config,repnum):
        """Initialize the Replica() object."""

        temp = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        self.repname = 'rep' + str(repnum).zfill(2)
        self.repdir = config.DATADIR + self.repname
        self.chain = Chain(config)
        self.mc = Monty(config,temp,self.chain)

