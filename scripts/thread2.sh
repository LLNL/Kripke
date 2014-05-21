#!/bin/bash

OUTNAME=thread2
NITER=10
ZONES=24,24,24
DLIST=2:8
GLIST=2:8

export OMP_NUM_THREADS=1
srun -n8 ./kripke --procs 2,2,2 --out ${OUTNAME}_1.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=2
srun -n4 ./kripke --procs 2,2,1 --out ${OUTNAME}_2.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=4
srun -n2 ./kripke --procs 2,1,1 --out ${OUTNAME}_4.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=8
srun -n1 ./kripke --procs 1,1,1 --out ${OUTNAME}_8.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

