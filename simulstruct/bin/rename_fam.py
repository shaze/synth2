#!/usr/bin/env python3

import sys
import random
import argparse

parser = argparse.ArgumentParser(description='Fiddle with the fam file')
parser.add_argument('input', metavar='str', type=str)
parser.add_argument('output', metavar='str', type=str)
parser.add_argument('--sex-error', dest='sex_error', action='store',default=0,type=float)

args = parser.parse_args()
f=open(args.input)
g=open(args.output,"w")
chrs="ABCDEFGHIJKLMNPQRSTVWXYZ"
chs=[c for c in chrs]
i=n=0
j=7
k=1
all_ids=set()
for line in f:
    data=line.split()
    while True:
        the_id="%s%s%s%s%02d"%(chrs[i],chrs[j],chrs[k],chrs[n%23],n)
        if the_id in all_ids:
           n=n+1
           j=(j+1)%23
        else:
           break
    all_ids.add(the_id)
    n=n+1
    if n%7==0: k=k+1
    if k>=len(chrs):k=0
    if k%13==0: j=j+1
    if j>=len(chrs):j=0    
    if j%3==0: i=i+1
    if i>=len(chrs):i=0
    if n%97==0:n=0
    data[0]=the_id
    data[1]=the_id
    if random.random()<args.sex_error and data[5] in ["1","2"]:
        data[5] = "0" if data[5]=="1" else "0"
    out="\t".join(data)+"\n"
    g.write(out)
g.close()
