#!/usr/bin/env python

import sys

class FastQfile(list):
    def __init__(self,fh):
        while fh:
            try:
                seq_id = fh.next().strip("\n")
                seq    = fh.next().strip("\n")
                plus   = fh.next().strip("\n")
                qual   = fh.next().strip("\n")
                assert seq_id.startswith('@')
                assert plus.startswith('+')
            except StopIteration:
                return
            except AssertionError:
                raise Exception("NOT A PROPERLY FORMATTED FASTQ FILE")
            except:
                raise
            else:
                self.append(FastQ(seq_id,seq,qual))
        return

class FastQ(object):
    def __init__(self,sqid,sq,ql):
        self.seqid = sqid
        self.seq = sq
        self.qual = ql
        return

    def __repr__(self):
        return "\n".join([self.seqid,self.seq,'+',self.qual])

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return "<FastQ entry of length %d>" % len(self)


if __name__=="__main__":
    if sys.argv[1] == "queryseq":
        query = sys.argv[2:]  if len(sys.argv) > 1 else None
        for s in FastQfile(sys.stdin):
            if query:
                printme = False
                for q in query:
                    if q in s.seq:
                        printme = True
                if printme:
                    print repr(s)
            else:
                print repr(s)

    elif sys.argv[1] == "queryid":
        with open(sys.argv[2]) as fh:
            ids = []
            for line in fh:
                ids.append(line.rstrip())

        for s in FastQfile(sys.stdin):
            if s.seqid.split()[0] in ids:
                print repr(s)
    else:
        print >> sys.stderr, "Unkown operation %s" % sys.argv[1]

