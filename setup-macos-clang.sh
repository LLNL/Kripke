#!/bin/bash

CXX=mpic++ CXXFLAGS="-Wall -O3 -g -std=c++14" cmake $1 -DENABLE_MPI=Off

