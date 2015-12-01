#!/bin/bash

module load cudatoolkit/7.5

cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain/chaos_5_x86_64_ib-ic16.cmake -DENABLE_CUDA=On -DENABLE_OPENMP=On $@


