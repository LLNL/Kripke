#!/bin/bash

module load cudatoolkit/7.5

cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain/chaos_5_x86_64_ib-ic16.cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=On $@


