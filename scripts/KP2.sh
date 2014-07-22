#!/bin/sh

NITER=10
ZONES=12,12,12
DLIST=1:12,2:6,3:4,4:3,6:2,12:1
GLIST=1:128,4:32,8:16,16:8,32:4,128:1
LEGENDRE=4
OUTFILE=KP2.out

./src/tools/kripke --procs 1,1,1 --out $OUTFILE --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST --legendre $LEGENDRE

