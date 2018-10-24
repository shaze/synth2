#!/usr/bin/env python3

import sys
import argparse
import os.path
import random
import bisect

def handleBatch(batchrange):
    blist = batchrange.split(",")
    batch = {}
    b=1
    for r in blist:
        m0,m1=list(map(int,r.split("-")))
        for i in range(m0,m1):
            batch[i]=b
        b=b+1
    return batch

def splitOnMissing(missfname,num_under,num_above,miss_cut):
    miss_cut = float(miss_cut)
    num_under = int(num_under)
    num_above = int(num_above)
    f = open(missfname)
    n=1
    misslist=[]
    f.readline()
    for line in f:
        line=line.strip().split()
        misslist.append((float(line[-1]),n))
        n=n+1
    misslist.sort()
    i=0
    for m, the_id in misslist:
        if m>miss_cut:
            where_miss=i
            break
        i=i+1
    if num_under > where_miss:
        sys.exit("Can't find %d inds with  missingness lt %f (%d)\n"%
                 (num_under,miss_cut,where_miss))
    the_batch = random.sample(misslist[:where_miss],num_under)+\
                random.sample(misslist[where_miss:],num_above)
    batches = { the_id:1 for (m, the_id) in the_batch }
    other = random.sample(list(range(n-1)),n//2)
    for i in other:
        miss, the_id = misslist[i]
        if the_id not in batches: batches[the_id]=2
    return batches

parser = argparse.ArgumentParser(description='Anonymise fam.')
parser.add_argument('--desex', dest='desex', action='store_true',default=False,
                    help='remove the sex')
parser.add_argument('--depheno', dest='depheno', action='store_true',default=False,
                    help='remove the pheno')
parser.add_argument('--batch', dest='batch', action='store',
                    help='list the batches')
parser.add_argument('--miss-batch', dest="missing", nargs=4, action='store', help='new fam file',default=[])
parser.add_argument('base', action='store', help='base of IDs')
parser.add_argument('orig_fam', action='store', help='original fam file')
parser.add_argument('new_fam', action='store', help='new fam file')
args = parser.parse_args()


if args.batch:
    batches = handleBatch(args.batch)
elif len(args.missing)==4:
    batches = splitOnMissing(*args.missing)

f = open(args.orig_fam)
g = open(args.new_fam,"w")
fname,suf= os.path.splitext(args.new_fam)

if args.batch or args.missing:
    h=open("%s.phe"%fname,"w")
    h.write("FID\tIID\tSEX\tBATCH\tCC\n")

ch1 = "JETJUNGLEVIP"
ch2 = "ABCEDEFGHIJKLMNOPQRSTUVWXYZ"

i=10000
n=1
for line in f:
    data = line.split()
    theid = args.base+ch1[i%len(ch1)]+ch2[i%len(ch2)]+ch1[(i+7)%len(ch1)]+"%03d"%((i%1000)//5)
    data[0]=theid
    data[1]=theid
    i += 1
    if args.depheno:
        data[-1]="0\n"
    if args.desex:
        data[-2]="0"
    g.write("\t".join(data)+"\n")
    if args.batch or args.missing:
        h.write("\t".join([data[0],data[1],data[4],str(batches.get(n,0)),data[5]])+"\n")
    n=n+1
g.close()
if args.batch or args.missing:
    h.close()
