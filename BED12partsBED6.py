#!/usr/bin/env python

__doc__ = '''extracts a specific part of a BED12 into BED6'''

import sys
from copy import deepcopy
# read arg
valid_args = {
    #'5prime': '5\' exon',
    #'3prime': '5\' exon',
    '5utr':'5\'-UTR (if available)',
    'exon': 'all exons',
    'cds': 'coding sequence',
    'intron':'5\'-UTR (if available)',
    '3utr':'3\'-UTR (if available)'
    ###'intergenic': 'intergenic sequence'
}

if len(sys.argv) < 2 or not set(sys.argv[1].split(',')).issubset(set(valid_args.keys())):
    print >> sys.stderr, '\n*** '+__doc__+' ***\n'
    print >> sys.stderr, "USAGE: %s <part>,<part>,... < IN > OUT\n" % (sys.argv[0])
    print >> sys.stderr, "Implemented are:"
    for k,v in valid_args.items():
        print >> sys.stderr, '\t',k,'\t',v
    sys.exit(1)
else:
    parts = sys.argv[1].split(',')


for line in sys.stdin:
    f = line.split()
    if len(f)<12:
        sys.exit('ERROR: only works with BED12 files')

    # get blocks
    rawBlocks = zip(map(int,f[10].rstrip(',').split(',')), map(int,f[11].rstrip(',').split(',')))
    blocks = {'5utr':[], 'cds':[], '3utr':[], 'exon':[], 'intron':[]}
    for i, b in enumerate(rawBlocks):
        blockStart = int(f[1])+b[1]
        blockEnd = blockStart+b[0]
        # strand define UTRs accoring to genes' strand
        left = '3utr' if f[5]=='-' else '5utr'
        rite = '3utr' if f[5]=='+' else '5utr'
        try:
            assert left != rite
        except AssertionError:
            print >> sys.stderr, 'ERROR: strand symbol in BED12 file not +/-'
            raise

        # introns
        if i>0:
            blocks['intron'].append((blocks['exon'][-1][1], blockStart))

        # exons
        blocks['exon'].append((blockStart,blockEnd))

        # split CDS and UTR
        if blockStart >= int(f[7]):
            blocks[rite].append((blockStart,blockEnd))
        elif blockEnd<=int(f[6]):
            blocks[left].append((blockStart,blockEnd))
        elif blockStart >= int(f[6]):
            if blockEnd <= int(f[7]):
                blocks['cds'].append((blockStart,blockEnd))
            else:
                blocks['cds'].append((blockStart,int(f[7])))
                blocks[rite].append((int(f[7]),blockEnd))
        else:  # blockStart<ATG
            if blockEnd<=int(f[7]):
                blocks[left].append((blockStart,int(f[6])))
                blocks['cds'].append((int(f[6]),blockEnd))
            else:
                blocks[left].append((blockStart,int(f[6])))
                blocks['cds'].append((int(f[6]),int(f[7])))
                blocks[rite].append((int(f[7]),blockEnd))
        # 2DO introns,exons,intergenic

    # print
    for p in parts:
        try:
            bl = blocks[p]
            assert len(bl) > 0
        except KeyError:
            raise
        except AssertionError:
            continue
        else:
            g = deepcopy(f)  # get a fresh copy each time (if multiple regions would be requested)

            # BED6
            bl.sort()
            for b in bl:
                print '\t'.join([f[0], str(b[0]), str(b[1]), f[3]+'.'+p, f[4], f[5]])

            '''
            # BED12 ### DOES NOT WORK WITH INTRONS and EXONS APPARENTLY
            g = deepcopy(f)  # get a fresh copy each time (if multiple regions would be requested)
            bl.sort()
            thinStart = min([b[0] for b in bl])
            thinEnd = max([b[1] for b in bl])
            g[1] = g[6] = str(thinStart)
            g[2] = g[7] = str(thinEnd)
            g[3] += '.'+p
            g[9] = str(len(bl))
            g[10] = ','.join(map(str,[x[1]-x[0] for x in bl]))
            g[11] = ','.join(map(str,[x[0]-thinStart for x in bl]))
            print '\t'.join(g)
            '''
