#!/usr/bin/env python3


import pandas as pd
import sys

if len(sys.argv)<=1:
    sys.argv=["removeRelInds.py","$missing","$ibd_genome","$outfname"]


EOL=chr(10)

imissf = pd.read_csv(sys.argv[1],delim_whitespace=True,index_col=["FID","IID"])
genomef = pd.read_csv(sys.argv[2],delim_whitespace=True,usecols=["FID1","IID1","FID2","IID2","PI_HAT"])

outf   =open(sys.argv[3],"w")


def getDegrees():
   elts = set(imissf.index.values)
   degd = {}
   rel  = {}
   for elt in elts: 
       degd[elt]=0
       rel[elt]=[]
   deg = pd.Series(degd)
   for i,row in genomef.iterrows():
       x=tuple(row[["FID1","IID1"]].tolist())
       y=tuple(row[["FID2","IID2"]].tolist())
       deg[x]=deg[x]+1
       deg[y]=deg[y]+1
       rel[x].append(y)
       rel[y].append(x)
   return rel, deg


rel, deg = getDegrees() 

remove=[]
candidates = deg[deg>=1].sort_values(ascending=False)
for i,c in candidates.iteritems():
    if deg[i]>0:
        remove.append(i)
        deg[i]=deg[i]-1
        for other in rel[i]:
            deg[other]=deg[other]-1

remove.sort()
outf.write(EOL.join(map (lambda x: "%s %s"%(x[0],x[1]),remove)))
outf.close()
        
