#!/bin/bash

module load cudatoolkit/7.5

. /usr/local/tools/dotkit/init.sh
use gcc-4.9.3p

cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain/chaos_5_x86_64_ib-gcc.cmake -DENABLE_CUDA=On -DENABLE_OPENMP=On $@


