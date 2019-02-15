#!/usr/bin/env python

__doc__ = '''
Replaces 3'UTR sequence in a BED12 file with another based on best overlap
'''
__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2014, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "0.2"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]

import sys
from optparse import OptionParser
import dcbio.parse.BEDfile

def ovp(a,b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def movp(a,b):
    '''overlapping function returns overhang (non-overlapping length)'''
    alen, blen = sum([abs(x[1]-x[0]) for x in a]), sum([abs(x[1]-x[0]) for x in b])
    overlap = 0
    spliceoverlap = 0
    for x in a:
        for y in b:
            overlap += ovp(x,y)
            if x[0] == y[0]:
                spliceoverlap += 1
            if x[1] == y[1]:
                spliceoverlap += 1
    return (alen-overlap, blen-overlap), (len(a)*2-spliceoverlap, len(b)*2-spliceoverlap)

def merge(iv):
    saved = list(iv[0])
    for st, en in sorted([sorted(t) for t in iv]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)


if __name__=="__main__":

    parser = OptionParser()
    parser.add_option("-1", "--primary", dest="primary", metavar="FILE", help="Primary annotation")
    parser.add_option("-2", "--secondary", dest="secondary", metavar="FILE", help="Secondary annotation (extends first)")
    parser.add_option("-o", "--overlaps", dest="overlaps", metavar="FILE", help="Overlapping IDs (eg from bedIntersect)")

    (options, args) = parser.parse_args()


    # read BED12 files
    with open(options.primary) as fh:
        one = BEDfile.BED(fh)

    # read second and index
    with open(options.secondary) as fh:
        two = BEDfile.BED(fh)
        two = two.getDict('name')

    # read overlaps
    overlaps = {}
    with open(options.overlaps) as fh:
        for line in fh:
            f = line.split()
            try:
                overlaps[f[0]].append(f[1])
            except:
                overlaps[f[0]] = [f[1]]

    # overlap each and decide on best overlap
    # attach UTR
    # report what has been merged
    for b in one:
        # split in blocks and mark stop codon
        stopcodon1 = int(b.thickEnd) if b.strand == '+' else int(b.thickStart)
        b.meta['blocks'] = b.blockSubsets()
        # statistics
        originalLength = len(b)
        sys.stderr.write('\t'.join(map(str,[b.chrom, b.chromStart, b.chromEnd, b.name])) + '\t')
        if b.strand == '+':
            sys.stderr.write(str(sum([x[1]-x[0] for x in b.meta['blocks'][2]])) + '\t')
        elif b.strand == '-':
            sys.stderr.write(str(sum([x[1]-x[0] for x in b.meta['blocks'][0]])) + '\t')
        else:
            raise Exception('Unkown strand')
        # get overlaps for all candidates (return mismatches, so 0 is the optimal case)
        stats = []
        try:
            overlaps[b.name]
        except KeyError:
            b.name += '|NO'
            print b
        else:
            # get overlap metrics
            for o in overlaps[b.name]:
                # check strand compatibility
                try:
                    assert two[o].strand == b.strand
                except:
                    raise Exception('Strands of overlapping features are not equal')

                two[o].meta['blocks'] = two[o].blockSubsets()
                # check if stop codon the same
                stopcodon2 = int(two[o].thickEnd) if two[o].strand == '+' else int(two[o].thickStart)
                stopdiff = abs(stopcodon1-stopcodon2)
                # check CDS and splice
                over, splice = movp(b.meta['blocks'][1],two[o].meta['blocks'][1])
                # aggregate stats for the 3'UTR (strand dependent)
                if b.strand == '-':
                    left_utrover, left_utrsplice = movp(b.meta['blocks'][0],two[o].meta['blocks'][0])
                    two[o].meta['stats'] = (stopdiff,splice,over,left_utrover,left_utrsplice)
                elif b.strand == '+':
                    rite_utrover, rite_utrsplice = movp(b.meta['blocks'][2],two[o].meta['blocks'][2])
                    two[o].meta['stats'] = (stopdiff,splice,over,rite_utrover,rite_utrsplice)
                else:
                    raise Exception('Unkown strand')
            # print results
            for o in sorted(overlaps[b.name], key=lambda x: two[x].meta['stats'][:3]):  # sort by first 3 statistics
                s = two[o].meta['stats']
                if s[0] == 0 and \
                    s[1][0] == 0 and \
                    s[2][0] == 0 and \
                    s[3][0] == 0 and \
                    s[4][1] > 0:  # check if 1. STOP codon 2. CDS splice/seq of a subset of refseq 3.  3'UTR is extended

                    # merge 3'UTR blocks
                    if b.strand == '-':
                        b.meta['blocks'][0] = two[o].meta['blocks'][0]
                    elif b.strand == '+':
                        b.meta['blocks'][2] = two[o].meta['blocks'][2]
                    else:
                        raise Exception('Unkown strand')

                    # rebuild everything from blocks
                    b.thickStart = min([ x[0] for x in b.meta['blocks'][1] ])
                    b.thickEnd = max([ x[1] for x in b.meta['blocks'][1] ])
                    b.chromStart = min([ x[0] for x in b.meta['blocks'][0] ] + [ b.thickStart ])
                    b.chromEnd = max([ b.thickEnd ] + [ x[1] for x in b.meta['blocks'][2] ])
                    offsets,lengths = [],[]
                    for exon in merge(b.meta['blocks'][0] + b.meta['blocks'][1] + b.meta['blocks'][2]):
                        offsets.append(exon[0]-b.chromStart)
                        lengths.append(exon[1]-exon[0])
                    b.blockCount = len(offsets)
                    b.blockSizes = ','.join(map(str,lengths))
                    b.blockStarts = ','.join(map(str,offsets))

                    # finalize
                    b.name += '|' + o
                else:
                    # no extension possible
                    b.name += '|NA'
                print b
                break  # only print first match-extensions
        if b.strand == '+':
            sys.stderr.write(str(sum([x[1]-x[0] for x in b.meta['blocks'][2]])) + '\t')
        elif b.strand == '-':
            sys.stderr.write(str(sum([x[1]-x[0] for x in b.meta['blocks'][0]])) + '\t')
        else:
            raise Exception('Unkown strand')
        sys.stderr.write(str(len(b)-originalLength) +'\n')

