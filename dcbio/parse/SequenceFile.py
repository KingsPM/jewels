#!/usr/bin/env python

import sys

class SequenceFile:
    def __init__(self,filename,filetype=None,parse=False):
        self.filename = filename
        self.filetype = filetype
        self.sequences = []
        if parse:
            self.parse()
        return

    def __str__(self):
        return '<SequenceFile %s>' % (self.filename)

    def __repr__(self):
        return self.file

    def _guess(self):
        with open(self.filename) as fh:
            for line in fh:
                if line.startswith('>'):
                    self.filetype = 'fasta'
                break  # only evaluate first line
            fh.close()
        return bool(self.filetype)

    '''internal parsers'''
    def _fasta(self,fi):
        n,s,c = [], '', []
        with open(self.filename) as fh:
            for line in fh:
                if line.startswith('>'):
                    if n and s:
                        yield Sequence(n,s,c)
                        n,s,c = [], '', []
                    # parse name (and coordinates)
                    f = line[1:].rstrip().split(' ')
                    n = f[0].split('|')
                    try:
                        c = f[1].split(':')
                    except:
                        c = []  #no coordinates
                elif line.startswith('#'):
                    pass  #comments
                else:
                    s += line.rstrip()
            fh.close()
        yield Sequence(n,s,c)

    '''public functions'''
    def basename(self):
        return '.'.join(self.filename.split('.')[:-1])

    def parse(self):
        if not bool(self.filetype):
            try:
                assert self._guess()
            except:
                print >> sys.stderr, "WARNING: Cannot guess file format of file %s" % (self.filename)
                sys.exit(1)
        # parse sequences (generator)
        self.sequences = [ s for s in getattr(self, '_' + self.filetype)(self.filename) ]
        return len(self.sequences)

    def write(self,overwrite=True):
        writeflag = 'w' if overwrite else 'a'
        if self.filetype.startswith('fa'):
            fh = open(self.filename,writeflag)
            for seq in self.sequences:
                print >> fh, repr(seq)
            fh.close()
        else:
            raise(BaseException('Unkown file format: ' + self.filetype))
        return


class Sequence:
    def __init__(self,n,s='',c=[],t={}):
        self.name = n
        self.sequence = s
        self.coordinates = c
        self.tags = t
        return

    def __len__(self):
        '''returns ungapped length'''
        return len(self.sequence) - self.sequence.count('-')

    def __str__(self):
        return '<Sequence %s %s> %s...' % ('|'.join(self.name), ':'.join(self.coordinates), self.sequence[:50])

    def __repr__(self):
        fastring = '>' + '|'.join(self.name)
        if self.coordinates:
            fastring += ' ' + ':'.join(self.coordinates)
        fastring += '\n'
        fastring += '\n'.join([ self.sequence[i:i+100] for i in range(0,len(self.sequence),100) ]) 
        return fastring


if __name__=="__main__":
    fi = sys.argv[1]
    print "FILE", fi
    sf = SequenceFile(fi)
    print "SequenceFile", sf
    sf.parse()
    print sf.sequences 
