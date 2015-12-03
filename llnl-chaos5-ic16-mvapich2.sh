#!/bin/bash

cmake -DCMAKE_TOOLCHAIN_FILE=cmake/Toolchain/chaos_5_x86_64_ib-ic16.cmake -DENABLE_OPENMP=On $@


