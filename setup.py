#!/usr/bin/env python

import os

project = 'kripke'

sys_type = 'unknown'
if 'SYS_TYPE' in os.environ:
	sys_type = os.environ['SYS_TYPE']


build_dir = "./%s-%s" % (project, sys_type)

print "SYS_TYPE:  %s" % sys_type
print "PROJECT:   %s" % project
print "BUILD DIR: %s" % build_dir

os.system("rm -rf %s" % build_dir)
os.makedirs(build_dir)
os.chdir(build_dir)
os.system("cmake ..")


