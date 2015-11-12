#!/bin/bash

. /usr/local/tools/dotkit/init.sh
use clang

CXX=mpiclang++ CXXFLAGS="-O3 -g -fopenmp -std=c++14 -DRAJA_PLATFORM_BGQ -DRAJA_COMPILER_CLANG -DKRIPKE_USE_OPENMP -stdlib=libc++" cmake $1 -DENABLE_MPI=Off

