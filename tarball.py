#!/usr/bin/env python

# This script generates the release tarball

import os

build_dir = "./kripke-tarball"

os.system("rm -rf %s" % build_dir)
os.makedirs(build_dir)
os.chdir(build_dir)
os.system("cmake ..")
os.system("make package_source")
os.system("mv *.tar.gz ..")
os.chdir("..")
os.system("rm -rf %s" % build_dir)
print "Tarball for Kripke generated"

