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
