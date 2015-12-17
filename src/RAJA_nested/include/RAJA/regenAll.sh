#!/bin/bash

./genView.py > View.hxx
./genLayout.py > Layout.hxx

./genForallN.py 2 > Forall2.hxx
./genForallN.py 3 > Forall3.hxx
./genForallN.py 4 > Forall4.hxx
