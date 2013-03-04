============================================================
 Parallelized Finite Difference Time Domain (FDTD) Solver
 Version 2.0, January 29 2011
 Author(s):  Michael Strickland
 Email:  strickland@fias.uni-frankfurt.de
============================================================

 DESCRIPTION 
------------------------------------------------------------ 

This code uses finite differences to solve the Schrodinger 
EQ in imaginary time for an arbitrary 3d potential.  It uses
the MPI (Message Passing Interface) standard.  The lattice
is equally partitioned into slices along the "x" direction.
Code can extract ground state and first few excited state 
wavefunction and energies.

REQUIREMENTS
------------------------------------------------------------

The MPI (Message Passing Interface) API must be be installed 
on your system. Currently tested against MPICH and OpenMPI. 
Can run on a single computational node or as many as you 
like. 

 COMPILING
------------------------------------------------------------

To compile, simply type "make" from the command line.

 USAGE
------------------------------------------------------------

All parameters are specified in the params file.  

To run:

make run1  # runs on 2 nodes - 1 master and 1 worker
make run   # runs on 3 nodes - 1 master and 2 workers
make run4  # runs on 5 nodes - 1 master and 4 workers

or generally

mpirun -np <1+Number of Worker Nodes> mpisolve

 DEBUGGING
------------------------------------------------------------

mpirun N -x DISPLAY run_gdb.csh mpisolve


 CONTRIBUTORS
-----------------------------------------------------------
Michael Strickland
Adrian Dumitru
Yun Guo

 LICENSE
------------------------------------------------------------

GNU General Public License (GPLv3)
See detailed text in license directory 

 ATTRIBUTION
------------------------------------------------------------

We ask that if you use this code for work which results in a
publication that you cite the following papers:

  M. Strickland and D. Yager-Elorriaga, "A Parallel 
    Algorithm for Solving the 3d Schrodinger Equation",
    Journal of Computational Physics 229, 6015 (2010).
    [http://arxiv.org/abs/0904.0939]

  A. Dumitru, Y. Guo, A. Mocsy and M. Strickland,
    "Quarkonium states in an anisotropic QCD plasma,"
    Phys. Rev. D 79, 054019 (2009).
    [http://arxiv.org/abs/0901.1998]

  A. Dumitru, Y. Guo, and M. Strickland, "The imaginary part
    of the static gluon propagator in an anisotropic 
    (viscous) QCD plasma, Phys. Rev. D 79, 114003, (2009).
    [http://arxiv.org/abs/0903.4703]

  M. Margotta, K. McCarty, C. McGahan, M. Strickland, 
    D. Yager-Elorriaga, Quarkonium states in a complex-
    valued potential, arXiv:1101.4651 (2011).
    [http://arxiv.org/abs/1101.4651]

