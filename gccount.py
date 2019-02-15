#!/usr/bin/env python

import sys

counts = { 'G':0, 'C':0, 'A':0, 'T':0, 'N':0 }
seqnum = 0

for line in sys.stdin:
    if line.startswith('>'):
        try:
            print >> sys.stderr, (counts['G']+counts['C']) / float(counts['A']+counts['T']+counts['G']+counts['C'])
        except:
            continue
        seqnum += 1
        print >> sys.stderr, seqnum, line,
        continue
    else:
        counts['G'] += line.count('G')
        counts['C'] += line.count('C')
        counts['A'] += line.count('A')
        counts['T'] += line.count('T')
        counts['N'] += line.count('N')

print (counts['G']+counts['C']) / float(counts['A']+counts['T']+counts['G']+counts['C'])
print counts
