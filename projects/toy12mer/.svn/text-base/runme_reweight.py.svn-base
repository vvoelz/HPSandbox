import os, sys, glob

if (1):
    sequences = ['KAAAAA']
    sequences.append('AKAAAA')
    sequences.append('AAKAAA')
    sequences.append('AAAKAA')
    sequences.append('AAAAKA')
    sequences.append('AAAAAK')

    sequences.append('FAAAAA')
    sequences.append('AFAAAA')
    sequences.append('AAFAAA')
    sequences.append('AAAFAA')
    sequences.append('AAAAFA')
    sequences.append('AAAAAF')

for seq in sequences:
    os.system('python reweightT.py microT.mtx microT_%s.mtx microstates.dat %s'%(seq,seq))
