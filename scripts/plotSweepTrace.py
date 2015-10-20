#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig = plt.figure()
ax = fig.add_subplot(1,1,1)


# Open trace files
for fname in sys.argv[1:]:
  # Extract the rank from the filename
  fparts = fname.split(".")
  fparts.reverse()
  rank = int(fparts[0])
  
  # Read the input file
  print "Reading data for rank %d" % rank
  with open(fname, "rb") as fh:
    for line in fh.readlines():
      line = line.rstrip()
      fields = line.split(" ")
      
      if fields[0] == 'sweep_kernel':
        color = 'blue'
      
      t0 = float(fields[1])
      t1 = float(fields[2])
      #print "%f - %f" % (t0, t1)
      
      ax.add_patch(patches.Rectangle( (t0, rank), t1-t0, 1.0, facecolor=color))    

plt.autoscale()

plt.show()

