#!/usr/bin/env python

__doc__='''merges sorted bed intervals and creates combined features
similar to mergeBed+cleanBEDnames but as below

mergeBed:
AAAAAAAAA
     BBBBBBBBB
CCCCCCCCCCCCCC (where C is A;B)

thisScript:
AAAAAAAAA
     BBBBBBBBB
AAAAACCCCBBBBB (where C is A;B)
'''

import sys
from collections import defaultdict

def cleanNames(s):
    return ';'.join(sorted(list(set(s.split(';')))))

class mBED(object):
    def __init__(self,fields):
        self.chr = fields[0]
        self.chrStart = int(fields[1])
        self.chrEnd = int(fields[2])
        self.segments = {}  # ends per name
        for n in list(set(fields[3].split(';'))):
            self.segments[n] = [ fields[0], int(fields[2]) ]  # chr, end (indexed by name)
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
            finishedSegments.append( [self.chr, self.chrStart, segmentTo, ";".join(sorted(self.segments.keys()))] )
            self.chrStart = segmentTo

            # delete ended segments to clean up
            ended = [ n  for n,s in self.segments.iteritems() if s[1] == segmentTo ]
            for e in ended:
                del self.segments[e]

        return finishedSegments


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
