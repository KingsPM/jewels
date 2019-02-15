#!/usr/bin/env python

import sys
import re
from hashlib import md5
from copy import deepcopy
from itertools import compress
from dcbio.misc.Term import Colors
from Bio import SeqIO
from collections import Counter

def triplets(s):
    if len(s) % 3 == 0:
        m0 = [m.start() for m in re.finditer('-', s[0::3])]
        m1 = [m.start() for m in re.finditer('-', s[1::3])]
        m2 = [m.start() for m in re.finditer('-', s[2::3])]
        if m0 == m1 == m2:
            return True
    return False

# alignment Parser (from filehandle or list of lines)
def nextMali(fh):
    # get function name from alntype
    lines = []
    alnNumber = None
    for line in fh:
        f = line.split()
        try:
            assert len(f) > 3 and not line.startswith('#')
        except:
            continue
        if alnNumber and f[0] != alnNumber:
            try:
                assert len(lines) > 0
            except:
                print >> sys.stderr, "ERROR: no data to build alignment"
                raise
            aln = Aln()
            aln.malis(lines)
            yield aln
            lines = []
        # append line to buffer
        lines.append(line)
        alnNumber = f[0]
    aln = Aln()
    aln.malis(lines)
    yield aln


# rudimentary seq class
class SimpleSeq:
    def __init__(self,inseq=None):
        self.seq = ''
        self.name = ''
        if inseq:
            self.seq = str(inseq.seq)
            self.name = str(inseq.name)
            self.id = str(inseq.id)
        return

    def __str__(self):
        return '>'+self.name+'\n'+'\n'.join([ self.seq[i::] for i in range(0,len(self.seq),70) ])+'\n' 

    def gapped(self):
        if '-' in self.seq:
            return True
        return False

    def removeGaps(self):
        print >> sys.stderr, self.seq
        self.seq = re.sub('-', '', self.seq)
        print >> sys.stderr, self.seq
        return

    def writefh(self,fh,seqformat='fasta',wrap=90):
        if seqformat.startswith('fa'):
            fh.write('>'+self.name+'\n')
            for i in range(0,len(self.name),wrap):
                fh.write(self.seq[i:wrap]+'\n')
        else:
            sys.exit('unknown file format => ' + seqformat + '\n')
        return

# simple alignment class
class Aln:
    def __init__(self):
        self.seqs = {}
        self.sorting = []
        self.inlength = None
        self.coloring = []
        self.flavour = '?'
        self.description = '?'
        return

    def truncateNames(self,to=8):
        # check if can be truncated
        this = self.sorting
        that = [ n[:to] for n in self.sorting ]
        '''replaces sequence names this with that'''
        try:
            assert len(that) == len(set(that))
        except:
            print >> sys.stderr, "## ERROR ## truncated names must be unambiguous and of equal length"
            raise
        # change sequence dictionary
        for i in range(len(self.sorting)):
            name = self.sorting[i]
            self.seqs[name[:8]] = self.seqs[name]
            del self.seqs[name]
            self.sorting[i] = name[:8]
        return

    def _flavouring(self):
        allseq = ''
        for s in self.seqs.values():
            allseq += s
        nucs = set(['A', 'T', 'G', 'C', 'N', '-'])
        if set(list(allseq.upper())).issubset(nucs):
            if triplets(allseq):
                self.flavour = 'c'
            else:
                self.flavour = 'n'
        else:
            self.flavour = 'P'
        return

    def _coloring(self):
        residues = [set(x) for x in zip(*list(self.seqs.values()))]
        for r in residues:
            p = ''
            if '-' in r:
                if len(r) == 1:
                    # allgap
                    pass
                elif len(r) == 2:
                    # same
                    p += Colors.GREEN
                else:
                    # variable
                    p += Colors.RED
            else:
                p += Colors.BOLD
                if len(r) == 1:
                    p += Colors.GREEN
                else:
                    # variable
                    p += Colors.RED
            self.coloring.append(p)
        return

    def containsSequence(self,seqname):
        for s in self.sorting:
            if s.startswith(seqname):
                return True
        return False

    def removeSequence(self,seqname):
        newsorting = []
        for s in self.sorting:
            if s.startswith(seqname):
                del self.seqs[s]
            else:
                newsorting.append(s)
        try:
            assert len(newsorting) != len(self.sorting)
        except:
            print >> sys.stderr, "### WARNING ### Sequence %s not in alignment %s" % (seqname, self.description)
            return
        else:
            self.sorting = newsorting
        return

    def removeGappedCols(self, colWidth=1):
        seqs = []
        for s in self.sorting:
            seqs.append(self.seqs[s])
        newseqs = [''] * len(seqs)
        self.inlength = 0
        colsToRemove = 0
        for p in zip(*seqs):
            if colsToRemove > 0:
                colsToRemove -= 1
            elif '-' in p:
                colsToRemove = colWidth - 1 # remove additional/followup columns if necessary
            else:
                self.inlength += 1
                for i in range(len(p)):
                    newseqs[i] += p[i]
        for i in range(len(newseqs)):
            self.seqs[self.sorting[i]] = newseqs[i]
        return

    def removeSoftmasking(self):
        self.seqs = { k: s.upper() for k, s in self.seqs.items() }
        return self

    def maskSTOP(self):
        stops = set(['TAG','TGA','TAA'])
        for s in self.seqs.keys():
            nsq = ''
            for i in range(0,len(self.seqs[s]),3):
                codon = self.seqs[s][i:i+3].upper()
                if codon in stops:
                    nsq += 'NNN'
                else:
                    nsq += codon
            self.seqs[s] = nsq
        return

    def _isGapless(self):
        for s in self.seqs.keys():
            if self.seqs[s].count('-') != 0:
                return False
        return True

    def windowstats(self, window=12, slide=3, seqs=None):
        '''returns identity statisitics for sligning windows'''
        # checks if gapless
        try:
            assert self._isGapless()
        except:
            print >> sys.stderr, "ERROR: there are gapped columns"
            raise
        # only calculate stats for sequence subset
        if seqs:
            try:
                statseqs = [ self.seqs[s] for s in seqs ]
            except:
                print >> sys.stderr, "### ERROR ### Sequence not found..."
                raise
        else:
            statseqs = self.seqs.values()

        # variable columns (ignoring N columns)_
        variable = [ 'N' in x and bool(float(len(set(x)))-1.9) or bool(len(set(x))-1) for x in zip(*statseqs) ]
        unmasked = [ 'N' not in x for x in zip(*statseqs) ]

        stats = []
        for i in range(0, len(variable)-window+1, slide):
            # mismatch count
            mismatchResidues = variable[i:i+window].count(True)
            unmaskedResidues = unmasked[i:i+window].count(True)
            # average GC content
            basecount = Counter(''.join([ x[i:i+window] for x in statseqs ]))
            GCcount = (basecount['G'] + basecount['C'])
            ATcount = (basecount['A'] + basecount['T'])
            avgGCcontent = GCcount / float(GCcount + ATcount)
            stats.append([i, window, mismatchResidues, unmaskedResidues, avgGCcontent])
        return stats

    def filter(self, threshold=6, window=12, slide=3, seqs=None):
        '''
        returns a filtered copy and the removed columns offsets
        only uses `seqs` to caluclate filter properties (by default uses all)
        '''
        # checks if gapless
        try:
            assert self._isGapless()
        except:
            print >> sys.stderr, "ERROR: there are gapped columns"
            raise
        # only calculate stats for sequence subset
        if seqs:
            try:
                statseqs = [ self.seqs[s] for s in seqs ]
            except:
                print >> sys.stderr, "### ERROR ### Sequence not found..."
                raise
        else:
            statseqs = self.seqs.values()
        # variable columns
        var = [ bool(len(set(x))-1) for x in zip(*statseqs) ]
        # create vector for comrpession
        varpos = []
        for i in range(0, len(var), slide):
            if var[i:i + window].count(True) > threshold:
                varpos += range(i, i+window)
        varpos = set(varpos)
        cvec = [ 0 if x in varpos else 1 for x in range(len(var)) ]
        assert len(cvec) == self.inlength
        # compress sequences (copy)
        cp = deepcopy(self)
        for s in cp.sorting:
            cp.seqs[s] = ''.join(compress(cp.seqs[s], cvec))
            cp.inlength = len(cp.seqs[s])
        return cp, [True if x==0 else False for x in cvec]

    def writefh(self, ofh, fmt="phy"):
        # calculate unique identifier string
        ident = md5(''.join(self.sorting)).hexdigest()
        if fmt.startswith('phy'):
            '''writes phylip format'''
            print >> ofh, "   " + str(len(self.seqs.keys())) + ' ' + str(self.inlength)
            for s in self.sorting:
                try:
                    self.seqs[s]
                except:
                    raise
                print >> ofh, s
                for i in range(0, len(self.seqs[s]), 60):
                    print >> ofh, self.seqs[s][i:i + 60]
        elif fmt.startswith('clustal'):
            '''writes clustalw format'''
            print >> ofh, "CLUSTAL W (1.82) multiple sequence alignment\n"
            constrack = self._conservation()
            namelen = max([len(s.split(':')[1]) for s in self.sorting])
            formatstring = '{0:<'+str(namelen)+'}      {1}'
            for i in range(0, self.inlength, 60):
                for s in self.sorting:
                    try:
                        self.seqs[s]
                    except:
                        raise
                    print >> ofh, formatstring.format(s.split(':')[1], self.seqs[s][i:i+60]) + " " + str(i+len(self.seqs[s][i:i+60]))
                print >> ofh, formatstring.format('', constrack[i:i+60]) + '\n'
        elif fmt.startswith('paml'):
            '''phylip format for paml with unique identifiers'''
            print >> ofh, "   " + str(len(self.seqs.keys())) + ' ' + str(self.inlength)
            idx = 0
            for s in self.sorting:
                try:
                    self.seqs[s]
                except:
                    raise
                idx += 1
                print >> ofh, ident[:6]+str(idx)
                for i in range(0, len(self.seqs[s]), 60):
                    print >> ofh, self.seqs[s][i:i + 60]
        elif fmt.startswith('mali'):
            for s in self.sorting:
                try:
                    self.seqs[s]
                except:
                    raise
                species, gene = s.split(':')
                print >> ofh, '\t'.join([self.description, species, gene, self.seqs[s]])
        elif fmt.startswith('fa'):
            for s in self.sorting:
                try:
                    self.seqs[s]
                except:
                    raise
                species, gene = s.split(':')
                print >> ofh, '>' + gene
                for w in range(0,len(self.seqs[s]),100):
                    print >> ofh, self.seqs[s][w:w+100]
        else:
            sys.exit("ERROR: Unkown output format specified")
        return ident

    def phylip(self, ifh):
        self.seqs = {}
        self.sorting = []
        self.inlength = None
        for line in ifh:
            if line.startswith(' '):
                pass
            elif "_" in line:
                current = line.rstrip()
                self.sorting.append(current)
                self.seqs[current] = ''
            else:
                self.seqs[current] += line.rstrip()
        self._checkLength()
        self._flavouring()
        assert self.flavour != 'c' or self.inlength % 3 == 0
        return

    def clustalw(self, ifh):
        self.seqs = {}
        self.sorting = []
        self.inlength = None
        for line in ifh:
            if len(line) < 2 or line.startswith(' ') or line.startswith('CLUSTAL'):
                pass
            elif "_" in line:
                current = line.rstrip()
                self.sorting.append(current)
            else:
                f = line.split()
                if len(f) == 2:
                    try:
                        self.seqs[f[0]] += f[1]
                    except:
                        self.sorting.append(f[0])
                        self.seqs[f[0]] = f[1]
                else:
                    print >> sys.stderr, f
        self._checkLength()
        self._flavouring()
        assert self.flavour != 'c' or self.inlength % 3 == 0
        return

    def fasta(self, ifh):
        self.seqs = {}
        self.sorting = []
        self.inlength = None
        for line in ifh:
            if line.startswith('>'):
                current = line[1:].rstrip()
                self.sorting.append(current)
                self.seqs[current] = ''
            else:
                self.seqs[current] += line.rstrip()
        self._checkLength()
        self._flavouring()
        assert self.flavour != 'c' or self.inlength % 3 == 0
        return

    def malis(self, ifh):
        # print >> sys.stderr, 'MALIS FROM', ifh
        self.seqs = {}
        self.sorting = []
        self.inlength = None
        for line in ifh:
            if line.startswith(' ') or len(line) < 4:
                pass
            else:
                f = line.split()
                self.description = f[0]
                current = ':'.join(f[1:3])
                self.sorting.append(current)
                self.seqs[current] = f[3]
        self._checkLength()
        self._flavouring()
        assert self.flavour != 'c' or self.inlength % 3 == 0
        return

    def maf(self, ifh):
        #s ci_aburtoni4.ab.gene.s187.37 0 207 + 207 MGFSQT
        self.seqs = {}
        self.sorting = []
        self.inlength = None
        for line in ifh:
            if line.startswith(' ') or line.startswith('#') or len(line) < 4:
                pass
            else:
                f = line.split()
                ff = f[1].split(':') if ":" in f[1] else f[1].split('.')
                current = ':'.join([ff[0], '.'.join(ff[1:])])
                self.description = ''
                self.sorting.append(current)
                self.seqs[current] = f[6]
        self._checkLength()
        self._flavouring()
        assert self.flavour != 'c' or self.inlength % 3 == 0
        return

    def _checkLength(self):
        self.inlength = len(self.seqs[self.sorting[0]])
        # check length and flush
        for s in self.seqs.keys():
            self.seqs[s] = self.seqs[s].replace(' ', '')
            try:
                assert not self.inlength or self.inlength == len(self.seqs[s])
            except:
                print >> sys.stderr, "ERROR: alignment length check failed"
                print >> sys.stderr, 'INLENGTH', self.inlength
                print >> sys.stderr, 'SEQ', s
                print >> sys.stderr, 'SEQLENGTH', len(self.seqs[s])
                seql = { k:len(x) for k,x in self.seqs.items() }
                for k,v in seql.items():
                    print >> sys.stderr, "\t", k, v
                raise
            else:
                self.inlength = len(self.seqs[s])
        # check indexing
        try:
            assert len(self.seqs.keys()) == len(self.sorting)
        except:
            sys.stderr.write('## FATAL ## sorting key error\n')
            raise
        return

    def _seqArray(self, seqs=None):
        '''reads/writes ordered array of sequences'''
        if seqs:
            # reapply sequences
            try:
                assert len(seqs) == len(self.sorting)
            except:
                sys.stderr.write('## FATAL ## cannot reapply array of sequences of different length\n')
                raise
            for i in range(len(seqs)):
                self.seqs[self.sorting[i]] = seqs[i]
            self._checkLength()
        # return sequence array
        newseqs = [self.seqs[s] for s in self.sorting]
        return newseqs

    def display(self, length=90, color=True, codons=True, blink=None, fh=sys.stderr):
        if color:
            self._coloring()
        l = {}
        colreset = [Colors.RESET] * length
        for i in range(0, self.inlength, length):
            for s in self.sorting:
                if i == 0:
                    l[s] = 0
                seq = self.seqs[s][i:i + length]
                beg = l[s]
                end = l[s] + len(seq) - seq.count('-')
                if color:
                    attribs = []
                    # color
                    attribs.append(self.coloring[i:i + length])
                    # codons
                    if codons and self.flavour == 'c':
                        codoncol = []
                        pos = beg
                        for c in range(len(seq)):
                            if seq[c] != '-':
                                if (pos / 3) % 2 == 0:
                                    codoncol.append(Colors.BGCYAN)
                                else:
                                    codoncol.append(Colors.BGBLACK)
                                pos += 1
                            else:
                                codoncol.append('')
                        attribs.append(codoncol)
                    if blink:
                        #blinking = [Colors.REVERSE if x else '' for x in blink[i:i+length]]
                        blinking = [Colors.UNDERSCORE + Colors.BLINK if x else '' for x in blink[i:i+length]]
                        attribs.append(blinking)
                    attribs.append(seq)
                    attribs.append(colreset)
                    seq = ''.join([''.join(x) for x in zip(*attribs)])
                    # print
                    fh.write(Colors.BLUE + Colors.BOLD + '{:<50}'.format(s) + Colors.RESET)
                    fh.write(Colors.CYAN + '{:>6}'.format(beg) + Colors.RESET)
                    fh.write(' {} '.format(seq))
                    fh.write(Colors.CYAN + '{:<6}'.format(end) + Colors.RESET)
                    fh.write('\n')
                else:
                    print >> fh,'{:<30}{:>6} {} {:<6}'.format(s, beg, seq, end)
                l[s] = end
            fh.write('\n')
        return

    def remove_empty_columns(self):
        # init
        oldseqs = self._seqArray()
        newseqs = []
        # check if removing gaps is needed and initialize
        for i in range(len(oldseqs)):
            newseqs.append('')
            if not '-' in oldseqs[i]:
                return
        # filter cols
        for i in range(len(oldseqs[0])):
            col = [s[i] for s in oldseqs]
            colchars = set(col)
            if len(colchars) != 1 or '-' not in colchars:
                for s in range(len(newseqs)):
                    newseqs[s] += oldseqs[s][i]

        # reapply sequences
        self._seqArray(newseqs)
        return self

    def _conservation(self):
        # init
        oldseqs = self._seqArray()
        constrack = ''
        # filter cols
        for i in range(len(oldseqs[0])):
            col = [s[i] for s in oldseqs]
            colchars = set(col)
            if len(colchars) == 1 and '-' not in colchars:
                constrack += '*'
            else:
                constrack += ' '
        return constrack

    def pairStats(self):
        '''returns alignment stats per pair'''
        results = []
        for i in range(len(self.sorting)):
            for j in range(i, len(self.sorting)):
                if i != j:
                    x = self.seqs[self.sorting[i]]
                    y = self.seqs[self.sorting[j]]
                    stats = [0, 0, 0, 0]  # residues, residues, overlapping, identity, splitstat, percentident
                    for p in range(self.inlength):
                        if x[p] != '-':
                            stats[0] += 1
                        if y[p] != '-':
                            stats[1] += 1
                            if x[p] != '-':
                                stats[2] += 1
                                if x[p] == y[p]:
                                    stats[3] += 1
                    # get split stat
                    stats.append(1.0-min(float(stats[2])/float(stats[0]),float(stats[2])/float(stats[1])))
                    stats.append(float(stats[3])/float(stats[2]) if stats[2] > 0 else 0.0)  #identity
                    # append results
                    results.append(tuple(stats))
        return results
