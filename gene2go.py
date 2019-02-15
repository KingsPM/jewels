#!/usr/bin/env python

import sys

try:
    fh = open(sys.argv[1],'r')
except:
    print >> sys.stderr, "Give input file please"

out = open(sys.argv[1] + '.gene2go','w')
for line in fh:
    f = line.split()
    if len(f) > 1:
        for c in f[1:]:
            g = c.split(',')
            for term in g:
                print >> out, f[0] + '\t' + term

out.close()
fh.close()


