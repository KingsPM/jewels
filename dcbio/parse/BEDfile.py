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
import md5

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


## BED parser
class BED(list):
    def __init__(self, fh, key=None):
        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.key = key or (lambda x: (x.chrom, x.chromStart, x.chromEnd))
        for line in fh:
            if line[0] == "#":
                continue
            if line.startswith('track'):
                continue
            self.append(BEDline(line))

        self.seqids = sorted(set(b.chrom for b in self))
        self.sort(key=self.key)
        return

    def get_order(self):
        return dict((f.name, (i, f)) for (i, f) in enumerate(self))

    def get_simple_bed(self):
        return [(b.chrom, i) for (i, b) in enumerate(self)]

    def getDict(self, keyattribute='name'):
        selfdict = {}
        for b in self:
            selfdict[getattr(b,keyattribute)] = b
        return selfdict

    # same as above but returns lists of isoforms
    def getDictList(self, keyattribute='name'):
        selfdict = {}
        for b in self:
            k = getattr(b,keyattribute)
            try:
                selfdict[k].append(b)
            except:
                selfdict[k] = [b]
        return selfdict

    def mergeByName(self, mergeby='name'):
        '''
        merge BED entries by name (create blocks)
        ignores blocks of merged entries
        '''
        raise Exception('Function not implemented yet...')
        #getattr(mergeby)
        return

# BED line
class BEDline(object):
    __slots__ = ("chrom", "chromStart", "chromEnd",
                 "name",  "score", "strand",
                 "thickStart", "thickEnd", "itemRgb",
                 "blockCount", "blockSizes", "blockStarts",
                 "fields","meta")


    def __init__(self, sline):
        args = sline.strip().split("\t")
        try:
            assert len(args) >= 3
        except:
            print >> sys.stderr, '##', sline, '##'
            raise Exception('BED line has less than 3 fields')
        self.fields = len(args)
        self.meta = {}
        self.chrom = args[0]
        self.chromStart = int(args[1])
        self.chromEnd = int(args[2])
        try:
            assert self.chromStart <= self.chromEnd
        except AssertionError:
            print >> sys.stderr, '## ERROR ##',sline
            raise Exception('End is less than Start coordinate')
        # optional fields (now with proper typecasting)
        if len(args) > 3:
            for i in range(3, len(args)):
                if i == 4:
                    # typecast
                    args[i] = float(args[i])
                if i in [6,7,9]:
                    # typecast
                    args[i] = int(args[i])
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

    def segments(self):
        '''returns all absoloute coordinate segments'''
        if self.fields >= 12:
            rawBlocks = zip(map(int,self.blockSizes.split(',')), map(int,self.blockStarts.split(',')))
            return [ (int(self.chromStart)+b[1], int(self.chromStart)+b[1]+b[0]) for b in rawBlocks ]
        return [ int(self.chromStart),int(self.chromEnd) ]


    def blockSubsets(self):
        '''returns absolute coordinate segments (thin, thick,thin)'''
        rawBlocks = self.segments()
        blocks = [[],[],[]] #leftThin,thick,riteThin
        for b in rawBlocks:
            blockStart = b[0]
            blockEnd = b[1]
            # split CDS and UTR
            if blockStart >= int(self.thickEnd):  #thinRite
                blocks[2].append((blockStart, blockEnd))
            elif blockEnd<=int(self.thickStart):  #thinLeft
                blocks[0].append((blockStart, blockEnd))
            elif blockStart >= int(self.thickStart):  #thick?
                if blockEnd <= int(self.thickEnd):  #thick
                    blocks[1].append((blockStart,blockEnd))
                else:  #thick,thinRite
                    blocks[1].append((blockStart,int(self.thickEnd)))
                    blocks[2].append((int(self.thickEnd),blockEnd))
            else:  #thinLeft,thick,?
                blocks[0].append((blockStart,int(self.thickStart)))
                if blockEnd<=int(self.thickEnd):  #thinLeft,thick
                    blocks[1].append((int(self.thickStart),blockEnd))
                else:  #thinLeft,thick,thinRite
                    blocks[1].append((int(self.thickStart),int(self.thickEnd)))
                    blocks[2].append((int(self.thickEnd),blockEnd))
        return blocks

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

    def _setDefaults(self):
        try:
            self.name
        except AttributeError:
            self.name = 'unknown'
        try:
            self.score
        except AttributeError:
            self.score = '1000'
        try:
            self.strand
        except AttributeError:
            self.strand = '+'
        try:
            self.thickStart
        except AttributeError:
            self.thickStart = self.chromStart
        try:
            self.thickEnd
        except AttributeError:
            self.thickEnd = self.chromEnd
        try:
            self.itemRgb
        except AttributeError:
            self.itemRgb = '255,0,0'
        if self.fields < 9:
            self.fields = 9
        return

    def setBlocks(self,blocks):
        # set blocks (makes BED9)
        self._setDefaults()
        # create block strings
        self.blockCount = str(len(blocks))
        blockSizes = []
        blockStarts = []
        for b in blocks:
            try:
                assert b[1] <= self.chromEnd and b[0] > self.chromStart
            except:
                raise Exception("Block is not start<end")
            blockSizes.append(b[1]-b[0])
            blockStarts.append(b[0]-self.chromStart)
        self.blockSizes = ",".join(map(str, blockSizes))
        self.blockStarts = ",".join(map(str, blockStarts))
        return

    def getBlocks(self, thick=False):
        '''returns blocks of BED line'''
        if self.fields >= 12:
            # BED12
            # starts = [ self.chromStart + s for s in map(int,self.blockStarts.split(',')) ]
            # length = map(int,self.blockSizes.split(','))
            # blocks = [ (starts[i],starts[i]+length[i]) for i in range(len(starts)) ]
            blocks = self.segments()
            if thick and (self.thickStart != self.chromStart or self.thickEnd != self.chromEnd):
                cutblocks = []
                ## cut to cds
                for b in blocks:
                    if b[1] < self.thickStart or b[0] > self.thickEnd:
                        continue
                    elif self.thickStart <= b[0] and b[1] <= self.thickEnd:
                        cutblocks.append(b)
                    else:
                        # intersecting block
                        left = b[0] if b[0] >= self.thickStart else self.thickStart
                        rite = b[1] if b[1] <= self.thickEnd else self.thickEnd
                        cutblocks.append((left,rite))


                return cutblocks
            else:
                return blocks
        else:
            if self.fields >= 7 and thick:
                return [ (self.thickStart, self.thickEnd) ]
            else:
                return [ (self.chromStart, self.chromEnd) ]

# BED line for fast overlapping
class mBED(object):
    def __init__(self,fields):
        self.chr = fields[0]
        self.chrStart = int(fields[1])
        self.chrEnd = int(fields[2])
        self.segments = {}  # ends per name
        # get score
        try:
            score = int(fields[4])
        except:
            score = 1
        # set identifier/name
        try:
            for n in list(set(fields[3].split(';'))):
                self.segments[n] = [ fields[0], int(fields[2]), score ]  # chr, end (indexed by name)
        except:
            n = md5.new(str(fields)).hexdigest()
            try:
                self.segments[n] = [ fields[0], int(fields[2]), score ]
            except:
                print >> sys.stderr, n
                print >> sys.stderr, fields
                print >> sys.stderr, score

                raise
        return

    def __str__(self):
        return '< '+'\t'.join([self.chr, str(self.chrStart), str(self.chrEnd), ','.join([n+':'+str(s[1]) for n, s in self.segments.iteritems()])]) + ' >'

    def _ovp(self, other):
        # overlap and bookended
        return True if self.chr == other.chr and other.chrStart <= self.chrEnd and other.chrEnd > self.chrStart else False

    def _add(self, other):
        # check if existing segments are extended (sanity check that they actually overlap!!!)
        for k, s in other.segments.iteritems():
            try:  # does segment exist?
                self.segments[k]
            except KeyError:  # create new segment
                self.segments[k] = s
            except:  # something else is wrong
                raise
            else:  # extend existing segment
                if s[1] > self.segments[k][1]:
                    self.segments[k][1] = s[1]
            finally:  # update chrEnd
                self.chrEnd = s[1] if s[1] > self.chrEnd else self.chrEnd
            return


    def _shift(self,toPosition):  #shift and prints to given position
        finishedSegments = []
        # sort segments by endpoints, print finished segements and deletes them
        while self.chrStart < toPosition:  # if there are multiple segments...

            # find lowest endpoint and collect names
            segmentTo = [ self.segments[k][1] for k in sorted(self.segments, key=self.segments.get) ]

            # set correct endpoint
            segmentTo = min(segmentTo) if min(segmentTo) <= toPosition else toPosition

            # make segment and update chrStart
            finishedSegments.append( [self.chr, self.chrStart, segmentTo, \
                ";".join(sorted(self.segments.keys())), sum([x[2] for x in self.segments.values()]) ] )
            self.chrStart = segmentTo

            # delete ended segments to clean up
            ended = [ n  for n,s in self.segments.iteritems() if s[1] == segmentTo ]
            for e in ended:
                del self.segments[e]

        return finishedSegments


if __name__ == "__main__":
    bed = BED(sys.stdin)
    for b in bed:
        print b

    #


