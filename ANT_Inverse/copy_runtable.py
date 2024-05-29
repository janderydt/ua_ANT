file_name = sys.argv[1]
  
with open(file_name, "r") as fi: # Open the file in read mode
  lines = fi.readlines()
  
with open(file_name, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("runtable="):
      runtable = ln.strip("runtable=").strip("\n")      

shutil.copyfile(runtable, runtable+".tmp")    