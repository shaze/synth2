#!/usr/bin/env python3

import sys

fname=sys.argv[1]
n=float(sys.argv[2])
orig=len(open(fname).readlines())
print(int(n*orig))
