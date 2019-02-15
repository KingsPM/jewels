#!/usr/bin/env python

import sys
import re
import copy

__doc__ = '''
MAF parser functions
'''
__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2012-2014, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "2.1"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


# CONVENIENCE FUNCTIONS
def findpos(s, f):
    return [m.start() for m in re.finditer(f, s)]

def splicepos(s, p):
    p.sort(reverse=True)
    for pos in p:
        s = s[:pos] + s[pos+1:]
    return s

def indexNot(l, v):
    for i in range(len(l)):
        if l[i] != v:
            return i
    return None

def _mergeSeqs(a, b):
    # -- -
    # X- X
    # XY ?
    # XN X
    try:
        assert len(a) == len(b)
    except:
        print >> sys.stderr, '\n', "A", len(a)
        print >> sys.stderr, "B", len(b)
        raise
    # remove softmasking
    out = ''
    for i in range(len(a)):
        if a[i] == '-' and b[i] != '-':
            out += b[i]
        elif a[i] != '-' and b[i] == '-':
            out += a[i]
        elif a[i] == '-' and b[i] == '-':
            out += '-'
        else:
            if a[i] !='N' and b[i] == 'N':
                out += a[i]
            elif a[i] =='N' and b[i] != 'N':
                out += b[i]
            elif a[i] != b[i]:
                out += '?'
            else:
                out += a[i]
    assert len(out) == len(a)
    return out

def _guessSeparator(l):
    for sep in [':','.']:
        f = l.split()
        sp = f[1].find(sep)
        if sp > 0:
            return sep
    raise Exception('Cannot guess species:sequence separator')

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

    def IDs(self,spec=True,chrom=False):
        '''returns IDs from block sequences'''
        ids = []
        for s in self.seqs:
            if spec:
                ids.append(s.species)
                if chrom:
                    ids[-1] += s._sep + s.seqid
            elif chrom:
                ids.append(s.seqid)
        return ids

    def pretty(self, wrap=100):
        '''return line wrapped'''
        printer = ''
        for i in range(0, len(self.seqs[0]), wrap):
            for s in self.seqs:
                printer += '{sq:<20} {pos:>9} {sequence}\n'.format(sq=s.id(20), pos=i, sequence=s.seq[i:i+wrap])
            printer += '\n'
        return printer

    def pearson(self,wrap=70):
        '''returns pearson/fasta'''
        printer = ''
        for s in self.seqs:
            printer += '>' + s.id() +'\n'
            for i in range(0, len(self.seqs[0]), wrap):
                printer += '{sequence}\n'.format(sequence=s.seq[i:i+wrap])
        return printer

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

    def remove_softmask(self):
        for s in self.seqs:
            s.seq = s.seq.upper()
        return

    def concat(self, other):
        '''concatenates sequences (in order, warns if different names)'''
        try:
            assert len(self.seqs) == len(other.seqs)
        except:
            print >> sys.stderr, "ERROR: cannot concatenate 2 alignments with different number of sequences"
            raise
        for i in range(len(self.seqs)):
            try:
                assert self.seqs[i].species == other.seqs[i].species
                assert self.seqs[i].seqid == other.seqs[i].seqid
            except:
                print >> sys.stderr, "ERROR: cannot merge alignment with different sequence order"
            self.seqs[i].seq += other.seqs[i].seq
            self.seqs[i].length += other.seqs[i].length
        return self

    def change_sequences(self, species=None, keep=False):
        try:
            assert isinstance(species, set)
        except:
            if species:
                species = set(species)
        newseq = []
        for s in self.seqs:
            if keep:
                if s.species in species:
                    newseq.append(s)
                else:
                    pass
            else:
                if s.species in species:
                    pass
                else:
                    newseq.append(s)
        self.seqs = newseq
        self._indexSeqs()
        return

    def keep_sequences(self, spec, chrom):
        '''variant of change_sequences, allowing for species and chromosome requirements'''
        try:
            assert isinstance(spec, set)
        except:
            spec = set(spec)
        try:
            assert isinstance(chrom, set)
        except:
            chrom = set(chrom)

        newseq = []
        for s in self.seqs:
            if s.species in spec and s.species + s._sep + s.seqid in chrom:
                newseq.append(s)
            else:
                pass
        self.seqs = newseq
        self._indexSeqs()
        return

    def remove_columns_with(self,rm):
        # init
        newseqs = []
        rmfound = False
        for i in range(len(self.seqs)):
            newseqs.append('')
            if rm in self.seqs[i].seq:
                rmfound = True
        if rmfound:
            # filter cols
            for i in range(len(self.seqs[0])):
                col = [s.seq[i] for s in self.seqs]
                colchars = set(col)
                if rm not in colchars:
                    for s in range(len(self.seqs)):
                        newseqs[s] += col[s]
            # reapply sequences
            for s in range(len(self.seqs)):
                self.seqs[s].seq = newseqs[s]
            return self

    def remove_empty_columns(self):
        # init
        newseqs = []
        for i in range(len(self.seqs)):
            newseqs.append('')
            if not '-' in self.seqs[i].seq:
                return
        # filter cols
        for i in range(len(self.seqs[0])):
            col = [s.seq[i] for s in self.seqs]
            colchars = set(col)
            if len(colchars) != 1 or '-' not in colchars:
                for s in range(len(self.seqs)):
                    newseqs[s] += col[s]
        # reapply sequences
        for s in range(len(self.seqs)):
            self.seqs[s].seq = newseqs[s]
        return self

    def informative_alignment(self):
        sys.stderr.write("## DEPRECTATION WARNING ## MAF.informative_alignment()\n")
        '''removes uninformative residues'''
        '''BE AWARE adjusts positions!'''
        blocks = self.blocks()
        if blocks:
            for i in range(1, len(blocks)):
                blocks[0].concat(blocks[i])
            return blocks[0]
        return None

    def count_changes(self, seq1, seq2):
        ###
        res = [ set(x) for x in zip(self.seqs[self.idx[seq1][0]].seq, self.seqs[self.idx[seq2][0]].seq) ]
        changed = 0
        for s in res:
            if len(s) > 1 and '.' not in s and '-' not in s:
                changed += 1
        return changed

    def fix_redundancy(self):
        if not self.ambiguous:
            return self
        # find species redundancy
        newseqs = []
        #newidx = {}
        for sp in self.idx.keys():
            if len(self.idx[sp]) > 1:
                # merge all to first
                for i in range(1, len(self.idx[sp])):
                    # merge with element 0
                    self.seqs[self.idx[sp][0]].seq = _mergeSeqs(self.seqs[self.idx[sp][0]].seq, self.seqs[self.idx[sp][i]].seq)
                    # change chromosome id
                    self.seqs[self.idx[sp][0]].seqid = 'scaffold_merge'
                    self.seqs[self.idx[sp][0]].length = len(self.seqs[self.idx[sp][0]].seq) - self.seqs[self.idx[sp][0]].seq.count('-')
                    self.seqs[self.idx[sp][0]].offset = 0
                    self.seqs[self.idx[sp][0]].seqsize = self.seqs[self.idx[sp][0]].length
            newseqs.append(self.seqs[self.idx[sp][0]])
            #newidx[sp] = [ len(newseqs)-1 ]
        self.seqs = newseqs
        #self.idx = newidx
        self._indexSeqs()
        return self

    def disambiguate_sequences(self, clean=None, strict=True):
        '''returns a dictionary of disambiguated sequences (strings)'''
        dseq = {}
        skip = set()
        for s in self.seqs:
            if s.species in dseq.keys():
                dseq[s.species] = _mergeSeqs(dseq[s.species],s.seq)
            else:
                dseq[s.species] = s.seq

        if clean:
            # remove empty columns in refspecies-
            if '-' in dseq[clean]:
                pos = findpos(dseq[clean],'-')
                for s in dseq.keys():
                    dseq[s] = splicepos(dseq[s],pos)

        if strict:
            #remove sequences with '?'
            for k in dseq.keys():
                if "?" in dseq[k]:
                    del dseq[k]
        return dseq

    def seqdic(self):
        '''return dictionary of sequences'''
        ret = {}
        for s in self.seqs:
            ret[s.id()] = s.seq
        return ret

    def consensus(self, names=None, consChar='.', spName='consensus', sqName='consensus', strip=False):
        # get sequence numbers
        cseqs = []
        if names:
            names = set(names)
            for n in names:
                cseqs += self.idx[n]
        else:
            cseqs = range(len(self.seqs))
        # create consensus
        seqs = [self.seqs[i] for i in cseqs]
        residues = [''.join(set(residue)) for residue in zip(*[x.seq for x in seqs])]
        newseq = ''
        for r in residues:
            if '-' in r:
                newseq += '-'
            elif len(r) > 1:
                newseq += consChar
            else:
                newseq += r
        newseqlen = len(newseq) - newseq.count('-')
        self.seqs.append(MAFseq('a\t' + spName + '.' + sqName + '\t0\t' + str(newseqlen) + '\t+\t' + str(newseqlen) + '\t' + newseq))
        self._indexSeqs()
        # strip sequences that were used for consensus
        if strip:
            self.change_sequences(names, keep=False)
        return self

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

    def speciesSequence(self,sep=':'):
        sp = set()
        for s in self.seqs:
            sp.add(sep.join([s.species,s.seqid]))
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
def next(fh,sep=':'):
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


# TEST
if __name__=="__main__":
    reqSpec = set(['Niloticus', 'Nyererei'])
    for block in next(sys.stdin):
        if set(block.species()).issuperset(reqSpec):
            #print "ORIGINAL"
            #print repr(block)
            block.change_sequences(reqSpec, keep=True)
            #print "STRIPPED"
            print repr(block)
            #l = len(block)
            block.remove_empty_columns()
            print "NO_EMPTY_COLS"
            print repr(block)
            #if len(block) != l:
            #    sys.exit(1)
        else:
            print reqSpec, block.species()

    sys.exit(5)

    for block in next(sys.stdin):
        print 'INTITAL'
        print repr(block)
        block.strip()
        print 'STRIPPED'
        print repr(block)
        bs = block.blocks()
        print 'SUB', len(bs)
        for b in bs:
            print repr(b)


