#!/bin/bash

#. /usr/local/tools/dotkit/init.sh
#use gcc-4.9.3p

CXX=mpiicpc-16.0.150 CXXFLAGS="-O3 -g -fopenmp -parallel-source-info=2 -unroll-aggressive -finline-functions -xavx -std=c++14 -DRAJA_PLATFORM_X86_SSE -DRAJA_COMPILER_ICC -DKRIPKE_USE_OPENMP" cmake $1 -DENABLE_MPI=Off

