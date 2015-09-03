AUTHORS
=======
Adam J. Kunen <kunen1@llnl.gov>
Peter N. Brown <brown42@llnl.gov>
Teresa S. Bailey <bailey42@llnl.gov>
Pete G. Maginot <maginot1@llnl.gov>


OVERVIEW
========
Kripke is a simple, scalable, 3D Sn deterministic particle transport code.  
Its primary purpose is to research how data layout, programming paradigms and 
architectures effect the implementation and performance of Sn transport.  
A main goal of Kripke is investigating how different data-layouts affect 
instruction, thread and task level parallelism, and what the implications are 
on overall solver performance.


REQUIRMENTS
===========
CMake 3.0 or later
C++ Compiler (g++, icpc, etc.)
MPI 1.0 or later


BUILDING AND RUNNING
====================
Kripke comes with a CMake based build system that requires out-of-source builds.


Quick Start
-----------

Step 1:  Create a build space
  (assuming you are starting in the Kripke root directory)
	mkdir build

Step 2: Run CMake in that build space
  cd kripke
	cmake ..

Step 3: Now make Kripke:
  make -j8
  
Step 4: Run the test suite to make sure it works
  make test
  
Step 5: Run Kripke's default problem:
  ./kripke
  



