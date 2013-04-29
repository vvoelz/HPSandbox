README: HPSandbox example programs
Last Updated:  March 2012

GETTING STARTED

Installation

Please define the HPSANDBOXHOME environment variable to be this directory 

In bash: $ export HPSANDBOXHOME=/Users/vince/scripts/HPSandbox/trunk/HPSandbox

COPYRIGHT
This python package is Copyright (C) 2007 Vincent Voelz <vvoelz@stanford.edu>.
Feel free to modify this code as needed, as long as you can keep it publicly available!

INTRO

HPSandbox is a set of Python objects that allow you to quickly write simple python scripts
to explore the two-dimensional HP lattice model of proteins of Chan and Dill.

PACKAGE CONTENTS

Chain.py           An object to represent the 2D HP lattice chain and its attributes, with method functions.
Config.py          A data structure to hold configuration parameters.
Monty.py           A collection of functions to perform Monte Carlo move-set operations on the Chain() object.
Replica.py         A container object, to hold the Chain() and Monty() objects
Trajectory.py      A set of functions for creating, reading, writing, and organizing trajectory files

/examples          A directory of example scripts
/sequences         Containing descriptions of the native states of foldable sequences:
                       /clist - contact state lists for chain lengths 10 through 21
                       /conf  - coordinates (conformations) for chain lengths 10 through 19    
                       COUNTS - text file counts of all unique (nonsymmetric) conformations for a given chain length

This package has been tested with Python 2.3 and 2.4.   Older/newer versions may work too, but haven't been tested.

         
SETUP

In order to get these example scripts to work correctly, you need to set up the following:

1) The HPSandbox directory (i.e python module) must be defined in your PYTHONPATH environment variable
2) In the /examples folder mcrex.conf file needs to be changed to reflect the absolute pathname
   of the sequences/clist/hp**  directory.



EXAMPLE SCRIPTS

Please see the /examples directory and the README file therein for some test scripts and examples showing
how to use the HPSandbox function.


DOCUMENTATION

The following documentation can be obtained using the pydoc standard module of python. For example:

>>> from HPSandbox import *
>>> import pydoc
>>> pydoc.help(Chain)


CLASSES
    Chain
    
    class Chain
     |  An object to represent the 2D HP lattice chain and its attributes, with method functions.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, config)
     |      Initialize the Chain object.
     |  
     |  contactstate(self)
     |      Return the contact state of the chain as a list of (res1,res2) contacts (duples),
     |      where the residue numbering starts at 0.
     |  
     |  grow(self)
     |      Add a new link onto the chain vector, updating the coords and viability correspondingly.
     |  
     |  hpstr2bin(self)
     |      Convert a string of type 'HPHPHPPPHHP' to a list of 1s and 0s.
     |  
     |  lastvec(self)
     |      Report the last entry on the list.
     |  
     |  nonsym(self)
     |      Many of the conformations are related by rotations and reflections.
     |      We define a "non-symmetric" conformation to have the first direction '0'
     |      and the first turn be a '1' (right turn)
     |      
     |      nonsym() returns 1 if the vec list is non-symmetric, 0 otherwise
     |  
     |  shift(self)
     |      Shifts the chain vector to the 'next' list, according to an enumeration scheme where
     |      the most distal chain vector is incremented 0->1->2->3.  After 3, the most distal vector
     |      element is removed, and the next most distal element is incremented.  If there are multiple
     |      "3" vectors, this process is done recursively.
     |      
     |      Example:
     |          [0,0,0,0] --> [0,0,0,1] 
     |          [0,0,1,2] --> [0,0,1,3] 
     |          [0,1,0,3] --> [0,1,1] 
     |          [0,3,3,3] --> [1]
     |      
     |      This operation is very useful for enumerating the full space of chain conformations.
     |      shift()  will also update the coords and the viability, accordingly.
     |      
     |      RETURN VALUES
     |      
     |          returns 1 if its the last possible "shift" --> i.e. if it's all 3's, the search is done
     |          returns 0 otherwise
     |  
     |  vec2coords(self, thisvec)
     |      Convert a list of chain vectors to a list of coordinates (duples).
     |  
     |  viability(self, thesecoords)
     |      Return 1 if the chain coordinates are self-avoiding, 0 if not.


    Config

    class Config
     |  A data structure to hold all the configuration data for an HP model calculation.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, filename=None)
     |      Initialize the configuration data structure, and default values.
     |  
     |  print_config(self)
     |  
     |  read_configfile(self, filename)
     |      Read in configuration parameters from file.  The file should have formatted rows
     |      consisting of two fields, separated by white-space (or any non-printing characters, like tabs):
     |      
     |      HPSTRING              PHPPHPPPHP 
     |      INITIALVEC           [0,0,0,0,0,0,0,0,0,0]
     |      ....


    DistRestraint

    class DistRestraint
     |  For now, this is a harmonic constraint over a squared distance D = d^2
     |  where D = sum_{i,j} d^2_ij over all contacts.
     |  
     |  Methods defined here:
     |  
     |  D(self, chain)
     |      Return the sum of squared-distances over the selected contacts.
     |  
     |  __init__(self, contacts, kspring)
     |      Initialize the DistRestraint object
     |  
     |  energy(self, chain)
     |      return the energy of the distance restraint


    Monty

    class Monty
     |  A collection of functions to perform Monte Carlo move-set operations on an HP lattice Chain object.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, config, temp, chain)
     |      Initialize the Monte Carlo object...
     |  
     |  energy(self, chain)
     |      Calculate potential energy of the chain.
     |  
     |  metropolis(self, replica)
     |      Accept Chain.nextvec over Chain.vec according to a Metropolis criterion.
     |  
     |  move1(self, replica)
     |      Apply moveset 'MC1' to the chain:
     |      (i)  three-bead flips
     |      (ii) end flips
     |      
     |      REFERENCE: Dill and Chan, 1994, 1996.
     |  
     |  move2(self, replica)
     |      Apply moveset MC2 to the chain:
     |      (i)   three-bead flips
     |      (ii)  end flips
     |      (iii) crankshaft moves
     |      (iv)  rigid rotations
     |      
     |      REFERENCE:  Dill and Chan, 1994, 1996
     |  
     |  
     |  move3(self, replica)
     |      Apply moveset 'MC3' to the chain.
     |      This is just a simple set to change the direction of a single chain link.
     |      Example:
     |          [0,0,0,0,0] --> [0,0,1,0,0]
     |      where {0,1,2,3}={n,e,s,w} direction
     |      
     |      About 5% viable moves are expected.
     |  
     |  move4(self, replica)
     |      Apply moveset 'MC4' to the chain:
     |      This is another vert simple moveset, to just change one angle in a rigid rotation
     |      Like 'MS3', this generates about 5% viable moves.


    Replica

    class Replica
     |  A container object, to hold the Chain() and Monty() objects
     |  
     |  Methods defined here:
     |  
     |  __init__(self, config, repnum)
     |      Initialize the Replica() object.


    Trajectory

    class Trajectory
     |  A set of functions for creating, reading, writing, and organizing trajcetory files
     |  
     |  Methods defined here:
     |  
     |  __init__(self, replicas, config)
     |      Initialize the trajectory object
     |  
     |  cleanup(self, replicas)
     |      Write any remaining points in the trajectory and energy buffers to file,
     |      and close any open file handles.
     |  
     |  dump_enequeue(self, replica)
     |      Dumps the queued energy values to the respective files and clears them for further use.
     |  
     |  dump_trjqueue(self, replica)
     |      Dump the queue to the the respective files and clear them for future use.
     |  
     |  mkdir(self, pathname)
     |      Automatically create directory if it doesn't exist
     |  
     |  queue_ene(self, replica)
     |      Queue an energy value to the buffer, for writing to file.
     |  
     |  queue_trj(self, replica)
     |      Queue a trajectory point to the buffer for writing to file
     |  
     |  write_eneheader(self, filename, replica)
     |      Write column headers for the energy file.




Vincent Voelz
vvoelz@stanford.edu
