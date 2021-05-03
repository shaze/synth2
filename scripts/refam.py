#!/usr/bin/env python3

import sys
import argparse
import os.path
import random
import bisect
import pandas as pd
import numpy as np
row_choice="ABCDEFGH"
ch1 = "JETJUNGLEVIP"
ch2 = "ABCEDEFGHIJKLMNOPQRSTUVWXYZ"


parser = argparse.ArgumentParser(description='Anonymise fam.')
parser.add_argument('--desex', dest='desex', action='store_true',default=False,
                    help='remove the sex')
parser.add_argument('--depheno', dest='depheno', action='store_true',default=False,
                    help='remove the pheno')
parser.add_argument('--batch', dest='batch', action='store',
                    help='list the batches')
parser.add_argument('--skew', dest='skew', action='store',default=[],
                    help='plates to be skewed')
parser.add_argument('--sex-error', dest='sexerr', type=float, action='store',default=0,
                    help='what proportion of mix-up')
parser.add_argument('--miss-batch', dest="missing", nargs=4, action='store', help='new fam file',default=[])
parser.add_argument('base', action='store', help='base of IDs')
parser.add_argument('orig_fam', action='store', help='original fam file')
parser.add_argument('new_fam', action='store', help='new fam file')
args = parser.parse_args()



def getPlates(batches):
    plate = {}
    for bnum in batches.values():
        curr_plate = "WP65"+"%04d"%(bnum*3)
        curr_elt = 0
        plate_num= 0
        for i in batches.keys():
            if batches[i] != bnum: continue
            if curr_elt==96:
                plate_num=plate_num+1
                curr_elt=0
            row = row_choice[curr_elt//12]
            col = curr_elt%12+1
            plate[i]=("%s%03d"%(curr_plate,plate_num),"%s%02d"%(row,col))
            curr_elt=curr_elt+1
    return plate

def handleBatch(batchrange):
    blist = batchrange.split(",")
    batch = {}
    b=1
    for bnum, r in enumerate(blist):
        m0,m1=list(map(int,r.split("-")))
        for i in range(m0,m1):
            batch[i]=bnum
        b=b+1
    return batch


def splitOnMissing(missfname,num_under,num_above,miss_cut):
    miss_cut = float(miss_cut)
    num_under = int(num_under)
    num_above = int(num_above)
    f = open(missfname)
    n=0
    misslist=[]
    f.readline()
    for line in f:
        line=line.strip().split()
        misslist.append((float(line[-1]),n))
        n=n+1
    misslist.sort()
    i=0
    where_miss = bisect.bisect(list(map(lambda x:x[0],misslist)),miss_cut)
    if num_under > where_miss:
        sys.exit("Can't find %d inds with  missingness lt %f (%d)\n"%
                 (num_under,miss_cut,where_miss))
    the_batch = random.sample(misslist[:where_miss],num_under)+\
                random.sample(misslist[where_miss:],num_above)
    batches = { the_id:1 for (m, the_id) in the_batch }
    other = list(set(range(n))-set(the_batch))
    for i in other:
        miss, the_id = misslist[i]
        if the_id not in batches: batches[the_id]=2
    return batches

def newID(row):
    i = row.name
    return args.base+"%s%s%05d"%(ch1[i//3982],ch2[i//1001],i%1001)

def skewPlates(all_info,skews):
    all_info.sort_values(by=["col","Well"],inplace=True)
    print(all_info[all_info["SamplePlate"]=="WP650006007"].head())
    for plate in skews.split(","):
        last = all_info.loc[all_info['SamplePlate']==plate,["FID","IID","sex","batch","CC"]].iloc[-1]
        curr=all_info.loc[all_info['SamplePlate']==plate,["FID","IID","sex","batch","CC"]].shift(1)
        curr.iloc[0]=last
        all_info.loc[all_info['SamplePlate']==plate,["FID","IID",'sex',"batch","CC"]]=curr.astype(dtype={'batch':np.int64,'CC':np.int64,'sex':np.int64})
    print("--------",all_info[all_info["FID"]=="RXJB00008"])



def mixUp(all_info,num_swaps):
    for i in range(num_swaps):
        x1=random.randint(0,N)
        if all_info.loc[x1]['SamplePlate'] in args.skew: continue
        print(x1,all_info.loc[x1,'sex'],all_info.loc[x1,'fid'],all_info.loc[x1,'FID'])
        all_info.loc[x1,'sex']=3-all_info.loc[x1,'sex']

fam_df = pd.read_csv(args.orig_fam,delim_whitespace=True,header=None,names=['fid','iid','pat','mat','sex','CC'],dtype={"CC":np.int64,"sex":np.int64})
N = len(fam_df)
if args.batch:
    batches = handleBatch(args.batch)
elif len(args.missing)==4:
    batches = splitOnMissing(*args.missing)
else:
    batches = {i:1 for i in range(N)}

plates = getPlates(batches)

pldict = pd.DataFrame.from_dict({'batch':batches,'plate_well':plates})
print("OOOOOO",fam_df[fam_df["fid"]=='AB1055'])

pldict['SamplePlate']=pldict['plate_well'].apply(lambda z:z[0])
pldict['Well']=pldict['plate_well'].apply(lambda z:z[1])
pldict['col'] =pldict['Well'].str.slice(start=1)
pldict['batch']=np.int64(pldict['batch'])
all_info = fam_df.join(pldict)
all_info['batch']=np.int64(all_info['batch'])
all_info['FID']=all_info.apply(newID,axis=1)
all_info['IID']=all_info.apply(newID,axis=1)
print("111--------",all_info[all_info["FID"]=="RXJB00008"])


orig_index=all_info.index
skewPlates(all_info,args.skew)
print(all_info[all_info["SamplePlate"]=="WP650006007"].head())
mixUp(all_info,int(N*args.sexerr))
all_info = all_info.reindex(index=orig_index)

all_info[['FID','IID','pat','mat','sex','CC']].to_csv(args.new_fam,index=False,header=False,sep='\t')

fname,suf= os.path.splitext(args.new_fam)

all_info[['FID','IID','sex','batch','SamplePlate','Well','CC']].to_csv("%s.phe"%fname,index=False,sep='\t')



if args.depheno:
    all_info['CC']=0
if args.desex:
    all_info['sex']=0


