#!/usr/bin/env python3

import sys

fname = sys.argv[1]
out     = sys.argv[2]
index={}
fam = open("%s.fam"%fname);
i=0
for line in fam:
    data = line.rstrip().split()
    person = "%s_%s"%(data[0],data[1]);
    index[person]=str(i)
    i=i+1

    
cf= open("%s.cluster1"%fname)
g=open(out,"w")
for line in cf:
    data = line.rstrip().split()
    data = list(map(lambda x:index[x],data[1:]))
    g.write(" ".join(data)+"\n")
g.close()
