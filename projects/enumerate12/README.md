# README


## Description 

 
Usage: `python enumerate.py <configfile>`

Try:  `$ python enumerate.py enumerate.conf`

This program will read in an HP chain specified in the configure file,
and perform a full enumeration of conformational space.

The problem tabulates:

    1) the density of states (in energies/contacts)

    2) the number density of unique contact states, i.e. disjoint collections
       of microscopic conformations all sharing a unique set of interresidue contacts. 

These values are printed as output.


## The configuration file `enumerate.conf`

```
HPSTRING                HPHPHPHPPHPH
INITIALVEC              [0,0,0,0,0,0,0,0,0,0,0]
randseed                345
eps                     -3.2
NREPLICAS               1
REPLICATEMPS            [300.0]
EXPDIR                  ./enumerate_data
PRINTEVERY              1000
TRJEVERY                1000
ENEEVERY                1000
NATIVEDIR               ../../sequences/clist/hp12
```
