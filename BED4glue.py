#!/usr/bin/env python

__doc__='''merges bookended BED entries with same name (input MUST be sorted to make this work)'''

import sys

last = None
for line in sys.stdin:
    f = line.split()
    if last and last[3] == f[3] and last[2] == f[1]:
        last[2:4] = f[2:4]
    else:
        if last:
            print '\t'.join(last)
        last = f
print '\t'.join(last)
