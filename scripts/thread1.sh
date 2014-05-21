#!/bin/bash

NITER=10
ZONES=12,12,12
DLIST=2:8
GLIST=2:8

export OMP_NUM_THREADS=1
srun -n1 ./kripke --procs 1,1,1 --out thread1_1.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=2
srun -n1 ./kripke --procs 1,1,1 --out thread1_2.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=4
srun -n1 ./kripke --procs 1,1,1 --out thread1_4.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=8
srun -n1 ./kripke --procs 1,1,1 --out thread1_8.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

