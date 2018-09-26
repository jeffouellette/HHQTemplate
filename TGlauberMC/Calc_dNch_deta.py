#!/usr/bin/python

import os, os.path, sys, string
from math import sqrt

outFile = open("dNch_deta.dat", "w")
outFile.write("#event  e0  dNch_dy\n")

startingEvent = 0
numEvents = 500

dNch_deta_averages = [0, 0, 0, 0, 0]
dNch_deta_errors = [0, 0, 0, 0, 0]
dNch_deta = [[], [], [], [], []]
num_events_processed = [0, 0, 0, 0, 0]
dNch_deta_max = 0.;

for eventNum in range(startingEvent, numEvents) :

  if not os.path.isfile("event"+str(eventNum)+"/logdir/b3d.log"): # check if b3d output exists
    continue

  #if not os.path.isfile("event"+str(eventNum)+"/logdir/b3d-analyze.log"): # check if b3d analyze output exists yet for this file (b3d-analyze.log exists iff b3d has finished processing)
    continue

  energyScalingCategory = int(eventNum/int(numEvents/5))
  energyScalingValue = 0.0016 + 0.0001*energyScalingCategory

  num_events_processed[energyScalingCategory] = num_events_processed[energyScalingCategory] + 1
  
  numlines = sum(1 for line in open("event"+str(eventNum)+"/logdir/b3d.log", "r") )
  file = open("event"+str(eventNum)+"/logdir/b3d.log", "r")
  fileEnum = enumerate(file)
  lineNum = None;
  for i,line in fileEnum:
    if "dNch/dy" in line:
      lineNum = i
  file.close()
  if lineNum is None:
    continue

  file = open("event"+str(eventNum)+"/logdir/b3d.log", "r")
  lines = file.readlines()
  line = lines[lineNum]

  dNch_dy_pos = line.find("dNch/dy") + 8;

  offset = 1
  dNch_dy = line[dNch_dy_pos:dNch_dy_pos+offset] 


  while dNch_dy[-1] != " ":
    offset = offset+1
    dNch_dy = line[dNch_dy_pos:dNch_dy_pos+offset]
  dNch_dy = dNch_dy[:-1] # drops the last whitespace character
  if dNch_deta_max < float(dNch_dy):
    dNch_deta_max = float(dNch_dy)

  print("event: " + str(eventNum) + ", e0= " + str(energyScalingValue) + ", dNch_dy= " + dNch_dy)
  dNch_deta[energyScalingCategory].append(float(dNch_dy))
  outFile.write(str(eventNum) + " " + str(energyScalingValue) + " " + dNch_dy + "\n") 
  file.close()
print("max dNch_dy = " + str(dNch_deta_max) + "\n")
outFile.close()


outFile2 = open("dNch_deta_summary.dat", "w")
outFile2.write("#e0 numEvents dNch_dy_avg dNch_dy_err\n")
for energyScalingCategory in range(5) :

  if num_events_processed[energyScalingCategory] != 0:
    energyScalingValue = 0.0016 + 0.0001*energyScalingCategory
    dNch_deta_averages[energyScalingCategory] = sum(dNch_deta[energyScalingCategory])/num_events_processed[energyScalingCategory]
    dNch_deta_errors[energyScalingCategory] = sqrt(sum([(dNch_deta_averages[energyScalingCategory] - x)**2 for x in dNch_deta[energyScalingCategory]]) / num_events_processed[energyScalingCategory]**2)
    outFile2.write(str(energyScalingValue) + " " + str(num_events_processed[energyScalingCategory]) + " " + str(dNch_deta_averages[energyScalingCategory]) + " " + str(dNch_deta_errors[energyScalingCategory]) + "\n")

print(dNch_deta_averages)
print(dNch_deta_errors)
outFile2.close()

