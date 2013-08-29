HPSandbox
===============

HPSandbox is a set of Python objects that allow you to quickly write simple python scripts
to explore the two-dimensional HP lattice model of proteins of Chan and Dill.


Last Updated:  March 2012

GETTING STARTED
===============

Installation
---------------

Please define the HPSANDBOXHOME environment variable to be this directory 

In bash: $ export HPSANDBOXHOME=/Users/vince/scripts/HPSandbox/trunk/HPSandbox

COPYRIGHT
---------------
This python package is Copyright (C) 2007 Vincent Voelz <vvoelz@stanford.edu>.
Feel free to modify this code as needed, as long as you can keep it publicly available!



PACKAGE CONTENTS
---------------

  * Chain.py           An object to represent the 2D HP lattice chain and its attributes, with method functions.
  * Config.py          A data structure to hold configuration parameters.
  * Monty.py           A collection of functions to perform Monte Carlo move-set operations on the Chain() object.
  * Replica.py         A container object, to hold the Chain() and Monty() objects
  * Trajectory.py      A set of functions for creating, reading, writing, and organizing trajectory files

/examples          A directory of example scripts
/sequences         Containing descriptions of the native states of foldable sequences:
                       /clist - contact state lists for chain lengths 10 through 21
                       /conf  - coordinates (conformations) for chain lengths 10 through 19    
                       COUNTS - text file counts of all unique (nonsymmetric) conformations for a given chain length

This package has been tested with Python 2.3 and 2.4.   Older/newer versions may work too, but haven't been tested.

         
SETUP

In order to get these example scripts to work correctly, you need to set up the following:

  *  The HPSandbox directory (i.e python module) must be defined in your PYTHONPATH environment variable
  *  In the /examples folder mcrex.conf file needs to be changed to reflect the absolute pathname
   of the sequences/clist/hp**  directory.



EXAMPLE SCRIPTS

Please see the /examples directory and the README file therein for some test scripts and examples showing
how to use the HPSandbox function.


DOCUMENTATION

The following documentation can be obtained using the pydoc standard module of python. For example:

 >>> from HPSandbox import *
 >>> import pydoc
 >>> pydoc.help(Chain)

Frequently Asked Questions (FAQ)
===============

What can HPSandbox do?
----------------------
HPSandbox can either 1) enumerate, or 2) perform Monte Carlo "dynamics" for 2-dimensional, square-lattice "bead-on-a-string" type chains.

How long a chain can I simulate?
----------------------
It depends on how long you are willing to wait. For instance, all conformations of 16-mers can be enumerated in a few minutes on a typical personal computer. Each increase in chain length adds a factor of about 2.7 to the calculation.

Can I use other potentials besides the HP model?
----------------------
Sure! But you'll have to put it in yourself. The code is designed for HP sequences, so if you want to study a model using beads of only two flavors, it is easy to just modify theMonty.energy()class function. More complicated models would require a more thorough, but straightforward reworking of the code.

What are the reference(s) for the HP 2D Model?
----------------------
Lau, K.F. and K.A. Dill. A Lattice Statistical Mechanics Model of the Conformational and Sequence Spaces of Proteins. Macromolecules 22: 3986-3997 (1989).
Dill, K.A., S. Bromberg, K. Yue, K.M. Fiebig, D.P. Yee, P.D. Thomas, and H.S. Chan. Principles of Protein Folding - A Perspective From Simple Exact Models. Protein Science 4: 561-602, 1995.
Lau and Dill (1989) is the first use of the HP model, while Dill et al. (1995) is a more comprehensive review.

What movesets are used for the Monte Carlo routines in HPSandbox?
----------------------
The movesets are described in Dill et al. Protein Science 4: 561-602, 1995:

<img src="http://dillgroup.stonybrook.edu/images/code-and-toys/hp-sandbox/movesets.png">


Vincent Voelz
vvoelz@stanford.edu
