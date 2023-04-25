# Monte Carlo sampling with Replica Exchange


This program will perform Monte Carlo sampling of 

```
Usage: python mcrex.py <configfile>

Try:  $ python mcrex.py mcrex.conf

This program will read in an HP chain and run parameters specified in the configure file,
and perform a replica exchange Monte Carlo simulation.

For the example "mcrex.conf", an 11-mer sequence is simulated, and the program ends when
the native conformation (contact state) is found.
A directory of results is output to directory ./mcrex_data


