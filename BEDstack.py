#!/usr/bin/env python

__doc__='''
BED parser for overlapping and stacking BED segments (cumulative scoring and naming)

eg.
AAAAAAAAA
     BBBBBBBBB
        CCC
1---12-23-31-1  sum of scores (1 by default)
AAAAADDDEFFBBB  where D=AB, E=ABC, F=DB
'''

import sys
from collections import defaultdict
import md5
import dcbio.parse.BEDfile as BEDfile

if __name__=="__main__":
    e = None # chromosome buffer
    for line in sys.stdin:
        fields = line.split()
        # check
        try:
            assert int(fields[1]) <= int(fields[2])
        except:
            print >> sys.stderr, fields
            raise Exception("coordinate error (start bigger than end)")
        # skip zero length
        if int(fields[1]) == int(fields[2]):
            continue
        # build extensible BED instance (mBED)
        f = mBED(fields)
        if not e:
            e = f
        else:
            # sorting check
            try:
                if e.chr == f.chr:
                    assert e.chrStart <= f.chrStart
            except AssertionError:
                raise Exception("input has to be a sorted bed file")
            # overlap and extend
            if e._ovp(f):
                # shift and print up to added position, add segments
                for s in e._shift(f.chrStart):
                    print "\t".join(map(str,s))
                e._add(f)
            else:
                # print everthing and set new
                for s in e._shift(e.chrEnd):
                    print "\t".join(map(str,s))
                e = f
    # do last one
    for s in e._shift(e.chrEnd):
        print "\t".join(map(str,s))

