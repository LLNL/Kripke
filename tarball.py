#!/usr/bin/env python

# This script generates the release tarball

import os
import sys

build_dir = "./kripke-tarball"

symref = os.popen("git symbolic-ref -q HEAD").read().strip(" \n\r")
symref_l = symref.split('/')
branch = 'none'
if len(symref_l) > 0:
  branch = symref_l[len(symref_l)-1]
print("Branch: %s" % branch)

os.system("rm -rf %s" % build_dir)
os.makedirs(build_dir)
os.chdir(build_dir)
os.system("cmake .. -DCPACK_SOURCE_PACKAGE_FILE_NAME=kripke-%s" % branch)
os.system("make package_source")
os.system("mv *.tar.gz ..")
os.chdir("..")
os.system("rm -rf %s" % build_dir)
print "Tarball for Kripke generated"

