
import os,re

listoffiles = []
for i in os.listdir(os.getcwd()):
    if os.path.isfile(i):
        f = open(i)
        for line in f:
            if line.startswith('#include "'):
                incl = line.split('"')[1]
                subsystem = incl.split("/")[0]
                package = incl.split("/")[1]
                listoffiles.append(subsystem+"/"+package)

setoffiles = set(listoffiles)
for e in sorted(setoffiles):
    print('<use name="'+e+'"/>') 
