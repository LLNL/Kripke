#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(1,1,1)


data = []
T = []

# Open trace files
for fname in sys.argv[1:]:

  # Extract the rank from the filename
  fparts = fname.split(".")
  fparts.reverse()
  rank = int(fparts[0])
  
  # Read the input file
  #print "Reading data for rank %d" % rank
  with open(fname, "rb") as fh:
    for line in fh.readlines():
      line = line.rstrip()
      fields = line.split(" ")
      
      if fields[0] == 'sweep_kernel':
        color = 'blue'
            
      t0 = float(fields[1])
      t1 = float(fields[2])
      T.append(t0)
      T.append(t1)
      data.append( (t0, 1) )
      data.append( (t1, -1) )

t_min = min(T)
t_max = max(T)
print "Total time: %f seconds" % (t_max - t_min)

# Sort data based on timestamp
data.sort(key=lambda tup: tup[0])

# Compute a curve that shows concurrency
ax_time = []
ay_concur = []
concur = 0
t_last = t_min
ave_concur = 0.0
for i in data:
  concur += i[1]
  ax_time.append(i[0])
  ay_concur.append(concur)
  dt = (i[0] - t_last)
  ave_concur += concur * dt / (t_max-t_min)
  t_last = i[0]

print "Max concurrency: %d" % max(ay_concur)
print "Ave concurrency: %f" % ave_concur    
             
plt.plot(ax_time, ay_concur)
#plt.hist(ay_concur)

plt.autoscale()
plt.show()

