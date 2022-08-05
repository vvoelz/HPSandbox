# README

## HPSandbox example: Monte Carlo Replica Exchange sampling

```
python mcrex.py

    usage:  mcrex.py <configfile>
    Try:    mcrex.py mcrex.conf

This program will read in an HP chain and run parameters specified in the config file (`*.conf`),
and perform a replica exchange Monte Carlo simulation.
```

For this example "mcrex.conf", an 11-mer sequence is simulated using replica-exchange Monte Carlo.
The program ends when
the native conformation (contact state) is found.  The set of native contacts must be known 
ahead of time, and read in from a file. 

This files can be found in `../../sequences/clist/`, which contains information for all foldable
N-mer sequences from N=10 through N=21.  These files are compressed (`.tar.gz`) in the github repo:

```
% ls ../../sequences/clist
hp10.tar.gz hp11.tar.gz hp13.tar.gz hp15.tar.gz hp17.tar.gz hp19.tar.gz hp21.tar.gz
hp11        hp12.tar.gz hp14.tar.gz hp16.tar.gz hp18.tar.gz hp20.tar.gz
```

..with the exception of file `../../sequences/clist/hp11/PHPPHPHPPHH.clist` to be used with this example.

## Output

A directory of results is output to directory ./mcrex_data
```


## EXPECTED OUTPUT

```
% python mcrex.py mcrex.conf 

#--------------Reading non-default Config.py file...--------------#
Configuration parameters from DEFAULT:

HPSTRING                       'PHPPHPHPPHH'
INITIALVEC                     [1, 0, 1, 2, 1, 2, 1, 2, 3, 3]
randseed                       345
eps                            -2.0
RESTRAINED_STATE               [(1, 4), (6, 9)]
KSPRING                        0.0
NREPLICAS                      8
REPLICATEMPS                   [275.0, 300.0, 325.0, 350.0, 400.0, 450.0, 500.0, 600.0]
MCSTEPS                        500000
SWAPEVERY                      500
SWAPMETHOD                     'randompair'
MOVESET                        'MS2'
EXPDIR                         './mcrex_results'
PRINTEVERY                     100
TRJEVERY                       100
ENEEVERY                       100
NATIVEDIR                      '../..//sequences/clist/hp11'
STOPATNATIVE                   1
NATIVE CLIST: [(1, 4), (1, 10), (4, 9), (6, 9)]
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
Trajectory REPFILES: [<_io.TextIOWrapper name='./mcrex_results/data/workspace/0/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/1/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/2/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/3/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/4/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/5/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/6/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/7/mc.trj' mode='w' encoding='UTF-8'>]
Traceback (most recent call last):
  File "mcrex.py", line 187, in <module>
    if string.strip(config.MOVESET) == 'MS1':
AttributeError: module 'string' has no attribute 'strip'
(base) vv@cst14330 mcrex % python mcrex.py mcrex.conf

#--------------Reading non-default Config.py file...--------------#
Configuration parameters from DEFAULT:

HPSTRING                       'PHPPHPHPPHH'
INITIALVEC                     [1, 0, 1, 2, 1, 2, 1, 2, 3, 3]
randseed                       345
eps                            -2.0
RESTRAINED_STATE               [(1, 4), (6, 9)]
KSPRING                        0.0
NREPLICAS                      8
REPLICATEMPS                   [275.0, 300.0, 325.0, 350.0, 400.0, 450.0, 500.0, 600.0]
MCSTEPS                        500000
SWAPEVERY                      500
SWAPMETHOD                     'randompair'
MOVESET                        'MS2'
EXPDIR                         './mcrex_results'
PRINTEVERY                     100
TRJEVERY                       100
ENEEVERY                       100
NATIVEDIR                      '../..//sequences/clist/hp11'
STOPATNATIVE                   1
NATIVE CLIST: [(1, 4), (1, 10), (4, 9), (6, 9)]
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
	Initializing Chain.py object...
	creating Monty.py object....
Trajectory REPFILES: [<_io.TextIOWrapper name='./mcrex_results/data/workspace/0/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/1/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/2/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/3/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/4/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/5/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/6/mc.trj' mode='w' encoding='UTF-8'>, <_io.TextIOWrapper name='./mcrex_results/data/workspace/7/mc.trj' mode='w' encoding='UTF-8'>]
38 production steps
replica      viablesteps  steps        MOVEaccept   viableswaps  swaps        SWAPaccept   
0            17           37           0.000        0            0            0.000        
1            19           37           0.000        0            0            0.000        
2            33           37           0.000        0            0            0.000        
3            31           37           0.000        0            0            0.000        
4            10           37           0.000        0            0            0.000        
5            36           37           0.000        0            0            0.000        
6            36           37           0.000        0            0            0.000        
7            36           37           0.000        0            0            0.000        
NATIVE CLIST: [(1, 4), (1, 10), (4, 9), (6, 9)]
replica  foundnative  contact state chainvec    
0        1            [(4, 9), (6, 9)] [1, 0, 0, 1, 1, 3, 2, 0, 1, 2]
1        1            [(1, 4)]     [2, 3, 1, 2, 3, 2, 3, 2, 3, 3]
2        1            [(1, 6), (1, 10), (6, 9)] [1, 2, 2, 1, 0, 0, 2, 0, 3, 1]
3        1            [(1, 4), (1, 6), (1, 10), (6, 9)] [0, 1, 2, 0, 3, 3, 0, 3, 2, 1]
4        1            [(1, 4), (1, 6)] [2, 1, 2, 3, 0, 1, 0, 1, 1, 1]
5        1            [(1, 4), (1, 6), (1, 10), (4, 9)] [0, 0, 2, 0, 2, 1, 3, 0, 0, 2]
6        1            [(1, 4), (1, 10), (4, 9), (6, 9)] [2, 1, 1, 3, 1, 2, 1, 3, 0, 3]
7        1            []           [1, 1, 1, 1, 2, 3, 0, 0, 3, 0]
```


