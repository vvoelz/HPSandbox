# README

## HPSandbox example: enumeration


```
python enumerate.py 

usage:  python enumerate.py <configfile>
Try:    puthon enumerate.py enumerate.conf

This program will read in an HP chain specified in the configure file,
and perform a full enumeration of conformational space.

The problem tablulates:
    1) the density of states (in energies/contacts)
    2) the number density of unique contact states, i.e. disjoint collections
       of microscopic conformations all sharing a unique set of interresidue contacts. 

These values are printed as output.


