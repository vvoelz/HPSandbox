import os, sys, glob

if (0):

    seqfile = '../seqdata.dat'
    fin = open(seqfile,'r')
    lines  = fin.readlines()
    fin.close()

    """reminder: data looks like this:
LYYYYY	   0.685	 357287.43417965	  61939.82758913	True
LYYYAY	   0.567	 164106.39817577	  32979.09622386	True
LYKYAY	   0.112	  77667.54783974	  14829.58358881	False
GYYYAY	   0.379	  10206.69795796	   5076.87412632	False
SYYYAY	   0.344	   9706.24854717	   4523.37210220	False
LYYQAY	   0.373	 137334.54134824	  21459.69614431	False
LYYYWY	   0.738	 573992.60766291	  91416.45612787	True
LYDYWY	   0.239	 302322.17593117	  32476.71646430	False
LYVYWY	   0.822	 669367.91996936	 134736.83195328	True
..."""

    total = 0
    for line in lines:

        fields = line.split()
        stability = float(fields[1])
        if stability > 0.4:
            seq = fields[0]
            os.system('python reweightT.py microT.mtx microT_%s.mtx microstates.dat %s'%(seq,seq))
            total += 1

        if total > 100:
            break


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
