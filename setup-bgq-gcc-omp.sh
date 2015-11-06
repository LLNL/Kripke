#!/bin/bash

. /usr/local/tools/dotkit/init.sh
use bggcc-4.7.2

CXX=mpig++-4.7.2 CXXFLAGS="-O3 -g -fopenmp -std=gnu++11 -DKRIPKE_USE_OPENMP" cmake $1 -DENABLE_MPI=Off

