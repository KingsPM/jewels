#!/usr/bin/env python

'''
substitutes chromosome names according to list
'''

import sys

def readlinks(fi):
    links = {}
    with open(fi) as fh:
        for line in fh:
            f = line.split()
            links[f[0]] = f[1]
    return links

if __name__=="__main__":
    skipped = {}
    li = readlinks(sys.argv[1])
    for line in sys.stdin:
        f = line.split('\t')
        try:
            sys.stdout.write('\t'.join([li[f[0]]] + f[1:]))
        except KeyError:
            try:
                skipped[f[0]] += 1
            except:
                skipped[f[0]] = 1
                print >> sys.stderr, "WARNING: unlinked chromosome %s (skipped)" % (f[0])

    print >> sys.stderr, "Skipped BED lines:"
    for k in sorted(skipped.keys()):
        print >> sys.stderr, '\t' + k + '\t' + str(skipped[k])
