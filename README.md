Authors
=======
  * Adam J. Kunen <kunen1@llnl.gov>
  * Peter N. Brown <brown42@llnl.gov>
  * Teresa S. Bailey <bailey42@llnl.gov>
  * Peter G. Maginot <maginot1@llnl.gov>


License
=======
See included file NOTICE


Overview
========
Kripke is a simple, scalable, 3D Sn deterministic particle transport code.  Its primary purpose is to research how data layout, programming paradigms and architectures effect the implementation and performance of Sn transport.  A main goal of Kripke is investigating how different data-layouts affect instruction, thread and task level parallelism, and what the implications are on overall solver performance.

Kripkie supports storage of angular fluxes (Psi) using all six striding orders (or "nestings") of Directions (D), Groups (G), and Zones (Z), and provides computational kernels specifically written for each of these nestings. Most Sn transport codes are designed around one of these nestings, which is an inflexibility that leads to software engineering compromises when porting to new architectures and programming paradigms.

Early research has found that the problem dimensions (zones, groups, directions, scattering order) and the scaling (number of threads and MPI tasks), can make a profound difference in the performance of each of these nestings. To our knowledge this is a capability unique to Kripke, and should provide key insight into how data-layout effects Sn solver performance. An asynchronous MPI-based parallel sweep algorithm is provided, which employs the concepts of Group Sets (GS) Zone Sets (ZS), and Direction Sets (DS), borrowed from the Texas A&M code PDT.

As we explore new architectures and programming paradigms with Kripke, we will be able to incorporate these findings and ideas into our larger codes. The main advantages of using Kripke for this exploration is that it's light-weight (ie. easily refactored and modified), and it gets us closer to the real question we want answered: "What is the best way to layout and implement an Sn code on a given architecture+programming-model?" instead of the more commonly asked question "What is the best way to map my existing Sn code to a given architecture+programming-model?".


Mini App or Proxy App?
----------------------
Kripke is a Mini-App since it has a very small code base consisting of 4178 lines of C++ code (generated using David A. Wheeler's SLOCCount v2.26).

Kripke is a Proxy-App since it is a proxy for the LLNL transport code ARDRA.


Analysis
--------
A major challenge of achieving high-performance in an Sn transport (or any physics) code is choosing a data-layout and a parallel decomposition that lends itself to the targeted architecture. Often the data-layout determines the most efficient nesting of loops in computational kernels, which then determines how well your inner-most-loop SIMDizes, how you add threading (pthreads, OpenMP, etc.), and the efficiency and design of your parallel algorithms. Therefore, each nesting produces different loop nesting orders, which provides substantially different performance characteristics. We want to explore how easily and efficiently these different nestings map to different architectures. In particular, we are interested in how we can achieve good parallel efficiency while also achieving efficient use of node resources (such as SIMD units, memory systems, and accelerators).

Parallel sweep algorithms can be explored with Kripke in multiple ways. The core MPI algorithm could be modified or rewritten to explore other approaches, domain overloading, or alternate programming models (such as Charm++). The effect of load-imbalance is an understudied aspect of Sn transport sweeps, and could easily be studied with Kripke by artificially adding more work (ie unknowns) to a subset of MPI tasks. Block-AMR could be added to Kripke, which would be a useful way to explore the cost-benefit analysis of adding AMR to an Sn code, and would be a way to further study load imbalances and AMR effects on sweeps.

The coupling of on-node sweep kernel, the parallel sweep algorithm, and the choices of decomposing the problem phase space into GS's, ZS's and DS's impact the performance of the overall sweep. The tradeoff between large and small "units of work" can be studied. Larger "units of work" provide more opportunity for on-node parallelism, while creating larger messages, less "sends", and less efficient parallel sweeps. Smaller "units of work" make for less efficient on-node kernels, but more efficient parallel sweeps. 

We can also study trading MPI tasks for threads, and the effects this has on our programming models and cache efficiency.

A simple timer infrastructure is provided that measure each compute kernels total time.


Physical Models
---------------
Kripke solves a simplification of the Discrete Ordinance and Diamond Difference discretized steady-state linear Boltzmann equation:
where (in this case) S is the identity matrix, and the 3D Diamond-Difference is used for inverting the "streaming and collision" operator (omega dot
grad + sigma). The kernels LTimes and LPlusTimes provide the action of the matrices L and L+, respectively.
One iteration of the MPI algorithm combined with the on-node sweep kernel solves one iteration of:


Implementation
--------------
The source code is divided into individual ".cpp" files, roughly one per class or function.
The file "src/Kripke/Kernel.h" defined the Kernel interface and defined a factory function for creating objects. The implemented kernels are in
"src/Kripke/Kernel/*.cpp". This is where most of the optimizations and programming model changes will occur.
The file "src/Kripke/Sweep_Solver.cpp" provides the MPI parallel sweep algorithm.
The rest of the files define the data structures that support these kernels and algorithms.
For a given problem, a single User_Data object is created, which contains everything: basic global parameters, quadrature set, MPI
communicator object (for Sweep_Solver), and a Grid_Data object.
The Grid_Data object contains L and L+ matrices, variables in moments space (phi), and a set of Group_Dir_Set objects for each <GS,DS> pair.
It is the Group_Dir_Set objects which store the angular fluxes(psi) and right hand side (rhs). In our implementation, we always ensure that all
directions contained in a single Group_Dir_Set object have the same sweeping order.



Inputs and Outputs
------------------
The Kripke build system produces a single executable "kripke". All of the parameters are supplied via command line options.
The number of GS and G (groups per groupset) are specified with "–grp <GS>:<G>". For example "–grp 16:2" specifies 16 groupsets and 2
groups per set, for a total of 32 groups.
Directions are specified almost the same way, with "–dir <DS>:<D>". However DS represents the number of direction sets PER octant.
Therefore "–dir 2:4" represents 2 direction sets per octant, 4 directions per set, and a total of 64 directions.
Nestings are selected with the "–nest DGZ" option, where any of the six DGZ, DZG, GDZ, GZD, ZDG, ZGD are allowed.
Number of Legendre moments ("–legendre") are selectable from 0 up to specify the number of moments in L and L+.
Number of zones (for the entire domain) are selected with "–zones X,Y,Z", number of MPI tasks are selected with "–procs X,Y,Z".
Number of iterations can be selected with "–niter N".
Search spaces can be specified by supplying lists to three arguments "–dir", "–grp", and "–nest". The product of these sets defines an entire
search space, all combinations of which will be run.
For example: "kripke --grp 1:4,2:2,4:1 --dir 16:1 --nest GDZ,ZDG" will run 6 different search points.


Building and Running
====================

Kripke comes with a simple CMake based build system.

Requirements
------------
*  CMake 3.0 or later
*  C++ Compiler (g++, icpc, etc.)
*  MPI 1.0 or later



Quick Start
-----------
The easiest way to get Kripke running, is to directly invoke CMake and take whatever system defaults you have for compilers and let CMake find MPI for you.

*  Step 1:  Create a build space

   > (assuming you are starting in the Kripke root directory)
   > 
	 > mkdir build

*  Step 2: Run CMake in that build space
   
   > cd kripke
	 > cmake ..

*  Step 3: Now make Kripke:
   
   > make -j8
  
*  Step 4: Run the test suite to make sure it works
   
   > make test
  
*  Step 5: Run Kripke's default problem:
   
   > ./kripke
  

Building with setup.py
----------------------
This needs to be written



Running Kripke
==============


Command Line Options
--------------------
Command line option help can also be viewed by running "./kripke --help"

### Problem Size Options:

* --groups <ngroups>     Number of energy groups
                         Default:  --groups 32

* --legendre <lorder>    Scattering Legendre Expansion Order (0, 1, ...)
                         Default:  --legendre 4

* --quad [<ndirs>|<polar>:<azim>]
                         Define the quadrature set to use
                         Either a fake S2 with <ndirs> points,
                         OR Gauss-Legendre with <polar> by <azim> points
                         Default:  --quad 96

* --zones <x,y,z>        Number of zones in x,y,z
                         Default:  --zones 16,16,16


### On-Node Options:

* --nest <NEST>          Loop nesting order (and data layout)
                         Available: DGZ,DZG,GDZ,GZD,ZDG,ZGD
                         Default:   --nest DGZ

###Parallel Decomposition Options:

* --layout <lout>        Layout of spatial subdomains over mpi ranks
                         0: Blocked: local zone sets are adjacent
                         1: Scattered: adjacent zone sets are distributed
                         Default: --layout 0

* --procs <npx,npy,npz>  Number of MPI ranks in each spatial dimension
                         Default:  --procs 1,1,1

* --dset <ds>            Number of direction-sets
                         Must be a factor of 8, and divide evenly the number
                         of quadrature points
                         Default:  --dset 8

* --gset <gs>            Number of energy group-sets
                         Must divide evenly the number energy groups
                         Default:  --gset 1

* --zset <zx>:<zy>:<zz>  Number of zone-sets in x:y:z
                         Default:  --zset 1:1:1

###Solver Options:

* --niter <NITER>        Number of solver iterations to run
                         Default:  --niter 10

* --pmethod <method>     Parallel solver method
                         sweep: Full up-wind sweep (wavefront algorithm)
                         bj: Block Jacobi
                         Default: --pmethod sweep


### Output and Testing Options:

* --out <OUTFILE>        Optional output file (default: none)

* --test                 Run Kernel Test instead of solve



Test Suite
----------



Future Plans
============
Some ideas for future study:
*  Block AMR.
*  More FLOP intensive spatial discretizations such as DFEM's.
*  Programming model abstractions


Retirement
==========
Retirement of this Mini-App should be considered when it is no longer a representative of state-of-the-art transport codes, or when it becomes too cumbersome to adapt to advanced architectures. Also, at the point of retirement it should be clear how to design its successor.


Release
=======
LLNL-CODE-658597
