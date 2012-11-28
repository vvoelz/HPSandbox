# Config.py     
#       
# HISTORY:
#
#     5/25/2007 - cleanup for hp module, added read_configfile when initialized
#     9/--/2005 - now reads in a configfile a la ZIPSEARCH
#     9/19/2005 - added self.STOPATNATIVE
#       

import string
import pickle

class Config:
  """A data structure to hold all the configuration data for an HP model calculation."""

  def __init__(self, filename=None):
    """Initialize the configuration data structure, and default values."""

    # DEFAULT VALUES
    self.HPSTRING = 'PHPPHPHPPHH'             # 11mer that zips, but with an implied contact
    self.INITIALVEC = [1,0,1,2,1,2,1,2,3,3]   # the initial starting conformation

    # Important constants, energy parms
    self.randseed = 345                     # random seed for random.Random()
    self.k =  0.001987                      # (kcal/K.mol) Boltzmann's constant 
    self.T = 300.0                          # reference temperataure units Kelvin (K)
    self.eps = -5.0                         # energetic strength of each contact (in units kT)
    self.epsilon = self.eps*self.k*self.T   # energetic strength of each contact -- in units Joules (J/mol)
                                            # kT = ~8.314 J/mol

    self.RESTRAINED_STATE = [(1,4),(6,9)]   # the contact state to be harmonically restrained
                                            # (in squared-distance d^2 = D)
    self.KSPRING = 0.0*self.k*self.T        # D spring constant (in kcal/mol per D)

    # Monte Carlo and Replica Exchange parameters
    self.NREPLICAS = 8                      # The number of replicas
    self.REPLICATEMPS = [275.0, 300.0, 325.0, 350.0, 400.0, 450.0, 500.0, 600.0]
				            # NOTE:  len(REPLICATEMPS) must equal NREPLICAS
    self.MCSTEPS = 500000                   # The number of Monte Carlo steps to perform
    self.SWAPEVERY = 500000                 # The frequency with which to attempt replica swaps 
    self.SWAPMETHOD = 'random pair'         # The method by which to select pairs of replicas for swaps
                                            #     options: 'random pair', and 'neighbors'
    self.MOVESET = 'MS2'                    # The type of Monte Carlo moveset
                                            #     options:   'MS1', 'MS2', 'MS3', and 'MS4'
                                            #     Please so Monty.py for the definitions of each move set

    self.PRINTEVERY = 1000                  # Frequency (in MC steps) to print status info to screen
    self.TRJEVERY = 1000                    # Frequency (in MC steps) to write trajectory data to file
    self.ENEEVERY = 1000                    # Frequenct (in MC steps) to write energy data to file 
    
    self.STOPATNATIVE = 1                   # 1 to stop the simulation as soon as the native is found, 0 if not.


    # Trajectory data directory pathnames
    self.EXPDIR = '/home6/voelzv/mc/exp/noswaps5e6'
    self.SETUPDIR = self.EXPDIR + '/setup'
    self.ANALDIR = self.EXPDIR + '/anal'
    self.DATADIR = self.EXPDIR + '/data'
   
    # Directory pathname where files of type 'HPPHPHPHPPHP.clist' are located
    # each belonging to a foldable sequence, and containing the native contact list 
    self.NATIVEDIR = '/home6/voelzv/zippers/native/lattice/clist/hp13'
    
    if filename != None:
        self.read_configfile( filename )
        self.filename = filename
    else:
        self.filename = 'DEFAULT'


  def read_configfile( self,filename ):
    """Read in configuration parameters from file.  The file should have formatted rows
    consisting of two fields, separated by white-space (or any non-printing characters, like tabs):

    HPSTRING              PHPPHPPPHP 
    INITIALVEC           [0,0,0,0,0,0,0,0,0,0]
    ....
""" 
    print '\n#--------------Reading non-default Config.py file...--------------#'
    
    fin = open(filename)
    lines = fin.readlines()
    
    i = 0
    for i in range(0,len(lines)):
    
      fields = string.splitfields(lines[i])
      if len(fields)>0:
	
	# Parse the variable assignments from the dummy Config file
	if fields[0] == 'HPSTRING':
	    self.HPSTRING = fields[1]
    
	if fields[0] == 'INITIALVEC':
	    self.INITIALVEC = eval(string.joinfields(fields[1:]))
    
	if fields[0] == 'randseed':
	    self.randseed = eval(fields[1])
    
	if fields[0] == 'eps':
	    self.eps = eval(fields[1])
    
	if fields[0] == 'RESTRAINED_STATE':
	   self.RESTRAINED_STATE = eval(string.joinfields(fields[1:]))
    
	if fields[0] == 'KSPRING':
	   self.KSPRING = eval(fields[1])
    
	if fields[0] == 'NREPLICAS':
	    self.NREPLICAS = eval(fields[1])
    
	if fields[0] == 'REPLICATEMPS':
	    self.REPLICATEMPS = eval(string.joinfields(fields[1:]))
    
	if fields[0] == 'MCSTEPS':
	    self.MCSTEPS = eval(fields[1])
    
	if fields[0] == 'SWAPEVERY':
	    self.SWAPEVERY = eval(fields[1])
    
	if fields[0] == 'SWAPMETHOD':
	    self.SWAPMETHOD = string.joinfields(fields[1:])
    
	if fields[0] == 'MOVESET':
	    self.MOVESET = string.joinfields(fields[1:])
    
	if fields[0] == 'EXPDIR':
	    self.EXPDIR = string.joinfields(fields[1:])
    
	if fields[0] == 'PRINTEVERY':
	    self.PRINTEVERY = eval(fields[1])
    
	if fields[0] == 'TRJEVERY':
	    self.TRJEVERY = eval(fields[1])
    
	if fields[0] == 'ENEEVERY':
	    self.ENEEVERY = eval(fields[1])

	if fields[0] == 'NATIVEDIR':
	    self.NATIVEDIR = fields[1]

	if fields[0] == 'STOPATNATIVE':
	    self.STOPATNATIVE = eval(fields[1])

    
    # end of line-reading loop
    
    self.SETUPDIR = self.EXPDIR + '/setup'
    self.ANALDIR = self.EXPDIR + '/anal'
    self.DATADIR = self.EXPDIR + '/data'
    
    self.epsilon = self.eps*self.k*self.T


  def print_config(self):

    print 'Configuration parameters from %s:'%self.filename
    print
    print '%-30s %s'%('HPSTRING',repr(self.HPSTRING))
    print '%-30s %s'%('INITIALVEC',repr(self.INITIALVEC))
    print '%-30s %s'%('randseed',repr(self.randseed))
    print '%-30s %s'%('eps',repr(self.eps))
    print '%-30s %s'%('RESTRAINED_STATE',repr(self.RESTRAINED_STATE))
    print '%-30s %s'%('KSPRING',repr(self.KSPRING))
    print '%-30s %s'%('NREPLICAS',repr(self.NREPLICAS))
    print '%-30s %s'%('REPLICATEMPS',repr(self.REPLICATEMPS))
    print '%-30s %s'%('MCSTEPS',repr(self.MCSTEPS))
    print '%-30s %s'%('SWAPEVERY',repr(self.SWAPEVERY))
    print '%-30s %s'%('SWAPMETHOD',repr(self.SWAPMETHOD))
    print '%-30s %s'%('MOVESET',repr(self.MOVESET))
    print '%-30s %s'%('EXPDIR',repr(self.EXPDIR))
    print '%-30s %s'%('PRINTEVERY',repr(self.PRINTEVERY))
    print '%-30s %s'%('TRJEVERY',repr(self.TRJEVERY))
    print '%-30s %s'%('ENEEVERY',repr(self.ENEEVERY))
    print '%-30s %s'%('NATIVEDIR',repr(self.NATIVEDIR))
    print '%-30s %s'%('STOPATNATIVE',repr(self.STOPATNATIVE))

