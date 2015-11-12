#!/bin/bash

CXX=mpic++ CXXFLAGS="-Wall -O3 -g -std=c++14 -DRAJA_PLATFORM_X86_SSE -DRAJA_COMPILER_CLANG" cmake $1 -DENABLE_MPI=Off

