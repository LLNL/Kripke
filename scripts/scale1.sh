#!/bin/bash

NITER=10
ZONES=12,12,12
DLIST=1:4,2:2,4:1
GLIST=1:4,2:2,4:1

srun -n1 ./kripke --procs 1,1,1 --out scale1_1.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST
srun -n2 ./kripke --procs 2,1,1 --out scale1_2.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST
srun -n4 ./kripke --procs 2,2,1 --out scale1_4.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST
srun -n8 ./kripke --procs 2,2,2 --out scale1_8.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST
srun -n16 ./kripke --procs 4,2,2 --out scale1_16.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST
