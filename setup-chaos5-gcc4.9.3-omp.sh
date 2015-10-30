#!/bin/bash

. /usr/local/tools/dotkit/init.sh
use gcc-4.9.3p

CXX=mpig++ CXXFLAGS="-O3 -g -fopenmp -mtune=native -std=c++14 -DKRIPKE_USE_OPENMP" cmake $1 -DENABLE_MPI=Off

