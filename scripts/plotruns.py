#!/bin/env python

import matplotlib as mpl
import numpy as np
from pylab import *

files = []
varnames = []
data = {}

normalize = False
reverse_order = False


# Read in and parse all of the data the user specified
for pair in sys.argv[1:]:
  (varname, fname) = pair.split(':')
  with open(fname) as fh:
    files.append(fname)
    varnames.append(varname)
    
    # read in file
    fdata = []
    idx = 0
    maxval = 1e100
    for line in fh:
      if not line.startswith("RUN:"):
        continue

      line = line.lstrip("RUN:")

      # cleanup
      line = line.strip()
      
      # split line into parts
      pairs = line.split(" ")
      
      # store each pair
      h = {}
      h['idx'] = idx
      idx += 1
      for p in pairs:
        if p != '':        
          (k, v) = p.split('=')
          h[k] = v
      
      # add data point for file
      fdata.append(h)
      num_points = idx
      maxval = min(maxval, float(h[varname]))
      
    # normalize if needed
    if normalize:
      for pt in fdata:
        pt[varname] = float(pt[varname])/maxval
      
    # add file data points to large hash
    data[fname] = fdata



# Sort the data
order = range(0,num_points)
varname = varnames[0]
def datasort(idx):
  minval = 1e100
  for fname in files:
    val = data[fname][idx][varname]
    minval = min([minval, float(val)])
  return minval

order = sorted(order, key=datasort)

if reverse_order:
  order.reverse()

# Generate a plot
fileidx = 0
for fname in files:
  fdata = data[fname]
  
  # pull out data
  X = []
  Y = []
  xval = 0
  for idx in order:
    pt = fdata[idx]        
    X.append(xval)
    Y.append(pt[varnames[fileidx]])
    if xval <= 10:
      print "%s: %s, dir=%s:%s,  grp=%s:%s,   Solve=%s, LTimes=%s, LPlusTimes=%s, Sweep=%s" % (fname, pt['nest'], pt['D'], pt['d'], pt['G'], pt['g'], pt['Solve'], pt['LTimes'], pt['LPlusTimes'], pt['Sweep'])
    xval += 1
  print ""
    
  # plot this files data
  plot(X, Y)
  
  fileidx += 1  
show()

    
