
#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
Simple BED file parser
'''

deprecation = '''
SimpleBED.py is

DD  EEE PPP RRR EEE  CC  A  TTT EEE DD
D D E   P P R R E   C   A A  T  E   D D
D D EE  PPP RRR EE  C   A A  T  EE  D D
D D E   P   RR  E   C   AAA  T  E   D D
DD  EEE P   R R EEE  CC A A  T  EEE DD
'''
raise Exception(deprecation)

import os.path as op
import itertools
import sys
import collections

## BED parser
class Bed(list):

    def __init__(self, filename, key=None):
        self.filename = filename
        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.key = key or (lambda x: (x.seqid, x.start, x.accn))
        for line in open(filename):
            if line[0] == "#": continue
            if line.startswith('track'): continue
            self.append(BedLine(line))

        self.seqids = sorted(set(b.seqid for b in self))
        self.sort(key=self.key)

    def get_order(self):
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

    def get_simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]

class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'stuff'.
    __slots__ = ("seqid", "start", "end", "accn", "stuff")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1])
        self.end = int(args[2])
        self.accn = args[3]
        self.stuff = args[4:] if len(args) > 4 else None

    def __str__(self):
        s = "\t".join(map(str, [getattr(self, attr) \
                    for attr in BedLine.__slots__[:-1]]))
        if self.stuff:
            s += "\t" + "\t".join(self.stuff)
        return s

    def __getitem__(self, key):
        return getattr(self, key)

