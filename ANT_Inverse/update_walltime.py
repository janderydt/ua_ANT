import datetime
import time
import sys

file_name = sys.argv[1]
seconds_expired = sys.argv[2]
  
with open(file_name, "r") as fi: # Open the file in read mode
  lines = fi.readlines()
  
with open(file_name, "r") as fi: # Open the file in read mode
  walltime = []
  for ln in fi:
    if ln.startswith("walltime="):
      walltime.append(ln[9:17])
      tmp = time.strptime(walltime[0],'%H:%M:%S')
      tmp_sec = datetime.timedelta(hours=tmp.tm_hour,minutes=tmp.tm_min,seconds=tmp.tm_sec).total_seconds()
      new_sec = tmp_sec - int(float(seconds_expired))
      new_time = time.strftime('%H:%M:%S', time.gmtime(new_sec))

kk = 0
with open(file_name, "r") as fi: # Open the file in read mode    
  for ln in fi: 
    if ln.startswith("walltime_remaining="):
      lines[kk] = "walltime_remaining="+new_time
    else:
      kk += 1
    
with open(file_name, "w") as fi: # Open the file in write mode
  fi.write("".join(lines)) # Write the modified content to the file
