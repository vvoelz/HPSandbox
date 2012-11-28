#! /usr/bin/python
#
# Replica.py    

from Config import *
import Chain
import Monty
from Trajectory import *

import random
import string
import math
import sys
import os

class Replica:
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self,config,repnum):
        """Initialize the Replica() object."""
	
	temp = config.REPLICATEMPS[repnum]
	self.repnum = repnum
	self.repname = 'rep' + string.zfill(str(repnum),2)
	self.repdir = config.DATADIR + self.repname
	self.chain = Chain.Chain(config)
	self.mc = Monty.Monty(config,temp,self.chain)    
    
