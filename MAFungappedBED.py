#!/usr/bin/env python

__doc__ = '''
returns ungapped blocks from MAF as BED file
'''

import sys
import re
import copy
from optparse import OptionParser

__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2014, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Prototype"

# CLASSES
class MAFseq:
    def __init__(self, l, sep=':'):
        f = l.split()
        self._sep = sep if sep else _guessSeparator(l)
        sp = f[1].find(sep)
        self.species = f[1][:sp]
        self.seqid = f[1][sp+1:]
        self.offset = int(f[2])
        self.length = int(f[3])
        self._strand = f[4]
        self.seqsize = int(f[5])
        self.seq = f[6]

        # sanity check
        try:
            assert len(self.seq) - self.seq.count('-') == self.length
        except:
            print 'TOT', len(self.seq)
            print 'GAP', self.seq.count('-')
            print 'SEQ', len(self.seq) - self.seq.count('-')
            print 'LEN', self.length
            print self
            raise
        return

    def __len__(self):
        return len(self.seq)

    def __str__(self, chars=20):
        p = ['']
        if len(self.seq) <= chars:
            p[0] += self.seq
            p[0] += (chars-len(self.seq)+3)*' '
        else:
            p[0] += self.seq[:chars/int(2)]
            p[0] += '...'
            p[0] += self.seq[-1*chars/int(2):]
        p.append('%9d' % self.start())
        p.append('%9d' % self.end())
        p.append('%2d' % self.strand())
        p.append(self.species)
        p.append(self.seqid)
        return '\t'.join(p)

    def __repr__(self):
        return '\t'.join(['s', self.species+self._sep+self.seqid, str(self.offset), str(self.length), self._strand, str(self.seqsize), self.seq])

    def BED(self):
        return '\t'.join(map(str,[self.seqid, self.start(), self.end(), self.id(), '.', self._strand]))


    def id(self, length=None):
        ident = self.species + self._sep + self.seqid
        if not length:
            return ident
        elif len(ident)<length:
            return ident + (" " * (length-len(ident)))
        return ident[:length-3] + ('.' * 3)

    def start(self):
        '''returns alignment start'''
        return self.end() - self.length

    def end(self):
        '''retuns alignment end'''
        if self.strand() < 0:
            return self.seqsize - self.offset
        return self.offset + self.length

    def strand(self):
        if self._strand == '+':
            return 1
        elif self._strand =='-':
            return -1
        return 0

class MAF:
    def __init__(self, seqs):
        self.seqs = seqs
        self.idx = {}
        self.ambiguous = False
        self._indexSeqs()
        return

    def __str__(self):
        p = ['<MAF Alignment>']
        for s in self.seqs:
            p.append(str(s))
        return '\n'.join(p)

    def __repr__(self):
        p = ['a']
        for s in self.seqs:
            p.append(repr(s))
        p.append('')
        return '\n'.join(p)

    def __len__(self):
        return len(self.seqs[0])

    def _indexSeqs(self):
        self.idx = {}
        for i in range(len(self.seqs)):
            try:
                assert self.seqs[i].species not in self.idx.keys()
            except:
                #print >> sys.stderr, "## WARNING ## multiple sequences for same species (index incomplete)"
                self.ambiguous = True
                self.idx[self.seqs[i].species].append(i)
            else:
                self.idx[self.seqs[i].species] = [i]
        return self.ambiguous

    def seq_by_name(self, name):
        try:
            assert name in self.idx.keys()
        except:
            print >> sys.stderr, name
            print >> sys.stderr, self.idx.keys()
            raise
        try:
            assert self.seqs[self.idx[name][0]].species == name
        except:
            print >> sys.stderr, self.idx
            print >> sys.stderr, self.seqs
            raise
        return self.seqs[self.idx[name][0]]

    def species(self):
        sp = set()
        for s in self.seqs:
            sp.add(s.species)
        return sp

    def _residueGaps(self):
        '''returns offsets of first non-informative residue with x sequences'''
        # zip sequences
        seqs = [s.seq for s in self.seqs]
        return [r.count('-') for r in zip(*seqs)]

    def split(self, at):
        '''returns new instance before split point, at column'''
        newseqs = []
        for s in self.seqs:
            # left (trim length)
            ss = copy.deepcopy(s)
            ss.length -= len(ss.seq[at:]) - ss.seq[at:].count('-')
            ss.seq = ss.seq[:at]
            # right (change offset and length)
            s.offset += len(s.seq[:at]) - s.seq[:at].count('-')
            s.length -= len(s.seq[:at]) - s.seq[:at].count('-')
            s.seq = s.seq[at:]
            newseqs.append(ss)
        return MAF(newseqs)

    def strip(self, reqseq=None):
        if reqseq is None:
            reqseq = len(self.seqs)
        '''strips ends that contain partly informative parts'''
        # get start and end of informative block
        rg = self._residueGaps()
        try:
            first = rg.index(0)
            last = len(rg) - rg[::-1].index(0)
        except ValueError:
            print self.pretty()
            raise
        if last <= first:
            # non informative
            return False
        # get
        for s in self.seqs:
            doffset = first - s.seq[:first].count('-')
            dlength = doffset + len(s.seq[last:]) - s.seq[last:].count('-')
            s.seq = s.seq[first:last]
            s.length -= dlength
            s.offset += doffset
        return True

    def blocks(self, reqseq=None):
        if reqseq is None:
            reqseq = len(self.seqs)
        '''returns MAF instances of informative blocks with reqseq sequences'''
        selfcopy = copy.deepcopy(self)
        if not selfcopy.strip(reqseq):
            return []
        blocks = []
        # determine offset where it should be split (beginning of noninformative block)
        rg = selfcopy._residueGaps()
        fr = indexNot(rg, 0)
        while fr:
            newaln = selfcopy.split(fr)
            blocks.append(newaln)
            selfcopy.strip(reqseq)
            rg = selfcopy._residueGaps()
            fr = indexNot(rg, 0)
        blocks.append(selfcopy)
        return blocks

# generator parsing function
def MAFnext(fh,sep=':'):
    seqs = []
    for line in fh:
        if line.startswith('#'):
            continue
        elif line.startswith('a') or len(line.rstrip()) == 0:
            if len(seqs) > 0:
                yield MAF(seqs)
                seqs = []
        elif line.startswith('s'):
            seqs.append(MAFseq(line,sep))

def indexNot(l, v):
    for i in range(len(l)):
        if l[i] != v:
            return i
    return None


if __name__=="__main__":
    # read requirements/includes
    usage = "usage: %prog -s species <MAF>"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", dest="species", metavar="STRING", help="Species for which the BED file is built")
    (options, args) = parser.parse_args()
    if not options.species:
        print >> sys.stderr, usage
        sys.exit(1)

    # select read handle
    try:
        readhandle = open(args[0],'r')
    except:
        if len(args)==0 or args[0]=='-':
            readhandle = sys.stdin
        else:
            sys.stderr.write('ERROR: Cannot open file ('+args[0]+')')
            raise

    for block in MAFnext(readhandle,':'):
        if not options.species in block.species():
            continue
        blocks = block.blocks()
        for b in blocks:
            #print repr(b)
            print b.seq_by_name(options.species).BED()

