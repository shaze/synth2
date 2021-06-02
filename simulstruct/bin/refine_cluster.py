#! /usr/bin/env python3


import numpy as np
import pandas as pd
import sys


df = pd.read_csv(sys.argv[1],header=None, names=['FID','IID','cluster'],
                 delim_whitespace=True)

min_size = int(sys.argv[2])

raw_clusters = df.groupby('cluster')



cluster_lens = sorted([(len(x), c, x) for c, x in raw_clusters ])

print([len(x) for c,x in raw_clusters])

def output(cluster_num, cluster):
    f=open("cluster_%03d"%cluster_num,"w")
    for x in cluster:
        f.write("%s %s\n"%(x,x))
    f.close()


curr_cluster = []
num_c = 0
for size, orig_c, cluster in cluster_lens:
    curr_cluster = curr_cluster + cluster['FID'].tolist()
    if len(curr_cluster)>min_size:
        print(len(curr_cluster))
        output(num_c, curr_cluster)
        curr_cluster = []
        num_c = num_c + 1

        
