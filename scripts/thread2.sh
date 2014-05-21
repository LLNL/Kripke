#!/bin/bash
#MSUB -V -j oe
#MSUB -l partition=rzmerl
#MSUB -l nodes=1
#MSUB -l walltime=08:00:00
#MSUB -A bdivp

#
# This runs on one multi-core node using 16 cores split between either MPI
# or OpenMP threads.
#
# Search space looks at 128 directions and 32 groups.
#

date

OUTNAME=thread2_rzmerl
NITER=10
ZONES=32,32,32
DLIST=1:16,2:8,4:4,8:2,16:1
GLIST=1:32,2:16,4:8,8:4,16:2,32:1

export OMP_NUM_THREADS=1
srun -n16 ./kripke --procs 4,2,2 --out ${OUTNAME}_1.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=2
srun -n8 ./kripke --procs 2,2,2 --out ${OUTNAME}_2.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=4
srun -n4 ./kripke --procs 2,2,1 --out ${OUTNAME}_4.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=8
srun -n2 ./kripke --procs 2,1,1 --out ${OUTNAME}_8.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

export OMP_NUM_THREADS=16
srun -n1 ./kripke --procs 1,1,1 --out ${OUTNAME}_16.out --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST

date

