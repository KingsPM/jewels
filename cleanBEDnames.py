#!/usr/bin/env python

'''
removes name redundancy in BED file after using mergeBed
'''

import sys

for line in sys.stdin:
    f = line.split()
    f[3] = ';'.join(list(set(f[3].split(';'))))
    print '\t'.join(f)

