#!/usr/bin/env python

import sys
from copy import deepcopy

keys = set()
fh = open(sys.argv[1])
for line in fh:
    keys.add(line.rstrip())
fh.close()
notfound = deepcopy(keys)
print >> sys.stderr, "Read %d terms to search" % (len(keys))

col = int(sys.argv[2])-1
print >> sys.stderr, "Will check in column %d" % (col)

cnt = [0,0,0]
for line in sys.stdin:
    f = line.split()
    try:
        if f[col] in keys:
            cnt[0] += 1
            notfound.remove(f[col])
            sys.stdout.write(line)
        else:
            print f[col]
            cnt[1] += 1
    except:
        cnt[2] += 1
        continue
    print >> sys.stderr, cnt, '\r',
print >> sys.stderr, ' matched', cnt[0]
print >> sys.stderr, ' nomatch', cnt[1]
print >> sys.stderr, '  failed', cnt[2]
print >> sys.stderr, 'notfound', len(notfound)
