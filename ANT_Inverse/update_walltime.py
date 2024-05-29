import datetime
import time
import sys
import math

file_name = sys.argv[1]
seconds_expired = sys.argv[2]
  
with open(file_name, "r") as fi: # Open the file in read mode
  lines = fi.readlines()
  
with open(file_name, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("walltime="):
      walltime = ln.strip("walltime=").strip("\n")
      tmp_sec = int(walltime.split(":")[0])*3600+int(walltime.split(":")[1])*60+int(walltime.split(":")[2])
      new_sec = tmp_sec - int(float(seconds_expired))
      new_hours = math.floor(new_sec/3600)
      new_minutes = math.floor((new_sec-new_hours*3600)/60)
      new_seconds = math.floor(new_sec-new_hours*3600-new_minutes*60)
      new_hours = str(new_hours).zfill(2)
      new_minutes = str(new_minutes).zfill(2)
      new_seconds = str(new_seconds).zfill(2)
      new_time = new_hours+':'+new_minutes+':'+new_seconds

kk = 0
with open(file_name, "r") as fi: # Open the file in read mode    
  for ln in fi: 
    if ln.startswith("walltime_remaining="):
      lines[kk] = "walltime_remaining="+new_time+"\n"
    else:
      kk += 1
    
with open(file_name, "w") as fi: # Open the file in write mode
  fi.write("".join(lines)) # Write the modified content to the file
