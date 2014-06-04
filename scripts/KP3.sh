#!/bin/sh

NITER=10
ZONES=12,12,12
DLIST=1:12,2:6,3:4,4:3,6:2,12:1
GLIST=1:64,2:32,4:16,8:8,16:4,32:2,64:1
LEGENDRE=9
OUTFILE=KP3.out

./src/tools/kripke --procs 1,1,1 --out $OUTFILE --zones $ZONES --niter $NITER --dir $DLIST --grp $GLIST --legendre $LEGENDRE

