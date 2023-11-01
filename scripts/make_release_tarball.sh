#!/usr/bin/env bash

###############################################################################
# Copyright (c) 2014-23, Lawrence Livermore National Security, LLC
# and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

TAR_CMD=`which tar`
VERSION=`git describe --tags`
SHORTHASH=`git rev-parse --short HEAD`

git archive --prefix=kripke-${VERSION}-${SHORTHASH}/ -o kripke-${VERSION}-${SHORTHASH}.tar HEAD 2> /dev/null

echo "Running git archive submodules..."

p=`pwd` && (echo .; git submodule foreach --recursive) | while read entering path; do
    temp="${path%\'}";
    temp="${temp#\'}";
    path=$temp;
    [ "$path" = "" ] && continue;
    (cd $path && git archive --prefix=kripke-${VERSION}-${SHORTHASH}/$path/ HEAD > $p/tmp.tar && ${TAR_CMD} --concatenate --file=$p/kripke-${VERSION}-${SHORTHASH}.tar $p/tmp.tar && rm $p/tmp.tar);
done

gzip kripke-${VERSION}-${SHORTHASH}.tar
