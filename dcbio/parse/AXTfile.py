#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
Complete BED12 class and parser (derived from SimpleBED which is now deprecated)

    chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 

The 9 additional optional BED fields are:

    name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
    strand - Defines the strand - either '+' or '-'.
    thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
    thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
    itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
    blockCount - The number of blocks (exons) in the BED line.
    blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
    blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
'''

import sys
import re
from itertools import count
import dcbio.parse.SimpleAlignment

# regex for axt parsing
regex_coord = re.compile('^\d+\s\w+\s\d+\s\d+\s\w+\s\d+\s\d+\s[+-]\s\d+$')
regex_nt = re.compile('^[AGCTN-]+$', re.IGNORECASE)

## Binary search for BEDline
def BEDbinsearch(key,lst,checkSort=False):
    # check if sorted
    if checkSort and not all(lst[i] <= lst[i+1] for i in xrange(len(lst)-1)):
        lst.sort()
    try: mid = lst[len(lst)/2]
    except: return []
    else:
        ov = []
        if key.overlaps(mid):
            ov.append(mid)
            for k in lst[len(lst)/2+1:]:
                if key.overlaps(k):
                    ov.append(k)
                else:
                    break
            for i in range(len(lst)/2-1,-1, -1):
                if key.overlaps(lst[i]):
                    ov.append(lst[i])
                else:
                    break
            return ov
        else:
            return BEDbinsearch(key,lst[:len(lst)/2]) + BEDbinsearch(key, lst[len(lst)/2+1:])



class AXTindex(object):
    """Indexed access to a axt using overlap queries, requires an index file"""
    def __init__( self, axt_filename, index_filename=None, keep_open=False, species1 = None, species2=None, species_to_lengths=None, support_ids=False ):
        if index_filename is None: index_filename = axt_filename + ".index"
        self.indexes = interval_index_file.Indexes( filename=index_filename )
        self.axt_filename = axt_filename
        # nota bene: (self.species1 = species1 or "species1") is incorrect if species1=""
        self.species1 = species1
        if (self.species1 == None): self.species1 = "species1"
        self.species2 = species2
        if (self.species2 == None): self.species2 = "species2"
        self.species_to_lengths = species_to_lengths
        self.support_ids        = support_ids            # for extra text at end of axt header lines
        if keep_open:
            self.f = open( axt_filename )
        else:
            self.f = None

    def get( self, src, start, end ):
        intersections = self.indexes.find( src, start, end )
        return itertools.imap( self.get_axt_at_offset, [ val for start, end, val in intersections ] )

    def get_axt_at_offset( self, offset ):
        if self.f:
            self.f.seek( offset )
            return read_next_axt( self.f, self.species1, self.species2, self.species_to_lengths, self.support_ids )
        else:
            f = open( self.axt_filename )
            try:
                f.seek( offset )
                return read_next_axt( f, self.species1, self.species2, self.species_to_lengths, self.support_ids )
            finally:
                f.close()

## AXT parser
class AXT(list):
    def __init__(self, fh, key=None):
        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.key = key or (lambda x: (x[0].targetChrom, x[0]targetStart, x[-1].targetEnd))
        self.append(AXTchain())

        lastBlock = -1
        lines = []

        for line in fh:
            if len(lines) == 0:
                if line.startswith('# identity'):
                    lines.append(line)
                else:
                    # leading comment line
                    pass
            elif len(lines) == 1 and line.startswith('# coverage'):
                lines.append(line)
            elif len(lines) == 2 and regex_coord.match(line) is not None:
                if lastBlock + 1 != int(line.split(' ')[0]):
                    # new chain (store old one)
                    self.append(AXTchain())
                lastBlock = int(line.split(' ')[0])
                lines.append(line)
            elif len(lines) == 3 and regex_nt.match(line) is not None:
                lines.append(line)
            elif len(lines) == 4 and regex_nt.match(line) is not None:
                lines.append(line)
                # parse and reset line buffer
                self[-1].add_aln(AXTaln(lines))
                lines = []
            else:
                raise Exception('AXT parsing error')

        # sort chains by target chromosome
        self.seqids = sorted(set(b.targetChrom for b in self))
        self.sort(key=self.key)
        ## indexing?

        return
    '''
    def get_order(self):
        return dict((f.name, (i, f)) for (i, f) in enumerate(self))

    def get_simple_bed(self):
        return [(b.chrom, i) for (i, b) in enumerate(self)]

    def getDict(self, keyattribute='name'):
        selfdict = {}
        for b in self:
            selfdict[getattr(b,keyattribute)] = b
        return selfdict
        '''


'''
AXTchain class
Counts Generated chains (gives thema unique id)
'''
class AXTchain(list):
    _ids = count(0)

    def __init__(self):
        self.id = self._ids.next()
        return

    def add_aln(self,aln):
        try:
            assert aln.num == len(self.aln) and aln.num == self[-1].num + 1
        except AssertionError:
            raise Exception("Alignment doesn't fit in chain")
        else:
            self.append(aln)
        return

    def targetRange(self):
        return (self[0].targetChrom, self[0].targetStart, self[-1].targetEnd)

    def queryRange(self):
        return (self[0].queryChrom, self[0].queryStart, self[-1].queryEnd)


'''
AXTalignment block
5 line block of an axt alignment
'''
class AXTaln(object):
    __slots__ = ("targetChrom", "targetStart", "targetEnd",
                 "queryChrom", "queryStart", "queryEnd",
                 "strand", "identity", "coverage",
                 "alnNumber", "meta")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        try:
            assert len(args) >= 3
        except:
            raise Exception('BED line has less than 3 fields')
        self.fields = len(args)
        self.meta = {}
        self.chrom = args[0]
        self.chromStart = int(args[1])
        self.chromEnd = int(args[2])
        self.aln = SimpleAlignment()
        try:
            assert self.chromStart <= self.chromEnd
        except AssertionError:
            print >> sys.stderr, '## ERROR ##',sline
            raise Exception('End is less than Start coordinate')
        # optional fields
        if len(args) > 3:
            for i in range(3, len(args)):
                if i in [10,11]:
                    # remove trailing commas in malformatted BED12 (thanks NCBI for those crappy RefSeq BEDs)
                    args[i] = args[i].rstrip(',')
                setattr(self, BEDline.__slots__[i], args[i])
        return

    def __str__(self):
        fields = []
        for i in range(self.fields):
            try:
                fields.append(getattr(self, BEDline.__slots__[i]))
            except:
                break
        return "\t".join(map(str, fields))

    def __getitem__(self, key):
        return getattr(self, key)

    def __lt__(self, other):
        return (self.chrom, self.chromStart, self.chromEnd) < (other.chrom, other.chromStart, other.chromEnd)

    def __len__(self):
        '''returns real length of intervals'''
        if self.fields >= 12:
            # returns bed12
             return sum(map(int,self.blockSizes.split(',')))
        else:
            return self.chromEnd - self.chromStart
        
    '''
    def overlaps(self,other):
        if (self.chrom != other.chrom) or (self.chromStart > other.chromEnd) or (self.chromEnd < other.chromStart):
            return False
        else:
            return True

    def overlapfrac(self,others,thick=False):
        '''return overlapping fraction allway'''
        coords = {}
        # own track
        for x in self.getBlocks(thick):
            try:
                coords[x[0]] += 1
            except:
                coords[x[0]] =  1
            try:
                coords[x[1]] -= 1
            except:
                coords[x[1]] = -1
        # others
        for other in others:
            for x in other.getBlocks(thick):
                try:
                    coords[x[0]] += 1
                except:
                    coords[x[0]] =  1
                try:
                    coords[x[1]] -= 1
                except:
                    coords[x[1]] = -1
        # evaluate
        #levels = []
        start, level, overlap = None, 0, []
        for k in sorted(coords.keys()):
            level += coords[k]
            if start is not None:
                assert level < len(others) + 1 
                overlap.append((start,k))
                start = None
            elif level == len(others) + 1:
                start = k
            #levels.append(level)
        #print >> sys.stderr, "L>", levels
        overlaplen = sum([ x[1]-x[0] for x in overlap ])
        return tuple([ float(overlaplen)/float(self.__len__()) ] + [ float(overlaplen)/float(len(o)) for o in others ])
    '''

if __name__ == "__main__":
    bed = BED(sys.stdin)
    for b in bed:
        print b
