#!/bin/bash

. /usr/local/tools/dotkit/init.sh
use clang-4.7.0

CXX=mpiclang++ CXXFLAGS="-O3 -g -fopenmp -mavx -std=c++14 -DRAJA_PLATFORM_X86_SSE -DRAJA_COMPILER_CLANG -DKRIPKE_USE_OPENMP" cmake $1 -DENABLE_MPI=Off

