#!/bin/bash

#. /usr/local/tools/dotkit/init.sh
#use gcc-4.9.3p

CXX=mpiicpc-16.0.109 CXXFLAGS="-O3 -g -fopenmp -unroll-aggressive -finline-functions -axAVX -std=c++14 -DKRIPKE_USE_OPENMP" cmake $1 -DENABLE_MPI=Off

