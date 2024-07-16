import sys
import shutil

config_file_name = sys.argv[1]
runtable_old_suffix = sys.argv[2]
runtable_new_suffix = sys.argv[3]
  
with open(config_file_name, "r") as fi: # Open the file in read mode
  lines = fi.readlines()
  
with open(config_file_name, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("runtable="):
      runtable = ln.strip("runtable=").strip("\n")      

shutil.copyfile(runtable+runtable_old_suffix, runtable+runtable_new_suffix)    
