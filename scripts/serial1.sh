#!/bin/sh

NITER=10
ZONES=12,12,12

echo "Running KP0"
./src/tools/kripke --procs 1,1,1 --out KP0.out --zones $ZONES --niter $NITER --dir 1:12 --grp 1:64 --legendre 4
echo "Running KP1"
./src/tools/kripke --procs 1,1,1 --out KP1.out --zones $ZONES --niter $NITER --dir 1:32 --grp 1:64 --legendre 4
echo "Running KP2"
./src/tools/kripke --procs 1,1,1 --out KP2.out --zones $ZONES --niter $NITER --dir 1:12 --grp 1:256 --legendre 4
echo "Running KP3"
./src/tools/kripke --procs 1,1,1 --out KP3.out --zones $ZONES --niter $NITER --dir 1:12 --grp 1:64 --legendre 9

