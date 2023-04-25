#!/usr/env python

import numpy as np

eps = -3.2  # units kT

dos = np.array([11460, 2850, 574, 130, 22, 1])
print dos
F = eps*np.arange(len(dos)) - np.log(dos)   # units kT

print '#number of contacts   energy (kT)  number of conformations  F(kT)  F(kcal/mol at T=300K)'
for i in range(len(F)):
    print '%4d %8.3f %16d %8.3f %8.3f'%(i, i*eps, dos[i], F[i], F[i]*0.5959)
print
print 'Total number of microstates:', dos.sum()

