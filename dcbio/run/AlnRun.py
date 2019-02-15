#!/usr/bin/env python

'''Class for alignment run and parse routines'''

import sys,math,StringIO
import itertools
import os
import subprocess
import tempfile
import time
from Bio import SeqRecord
import bx.align.lav
import bx.align.maf

from dcbio.parse.SimpleAlignment import SimpleSeq, Aln
from dcbio.misc.TmpFiles import TmpFiles
import dcbio.parse.MAF



class TmpFiles:
    def __init__(self,number=0):
        self.handles = []
        for i in range(number):
            self._tmpfile()
        return

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for fh in self.handles:
            try:
                fh.close()
            except:
                continue
            try:
                os.unlink(fh.name)
            except:
                continue
        return

    def _tmpfile(self):
        self.handles.append(tempfile.NamedTemporaryFile(delete=False))
        return self.handles[-1]


class Align:
    def __init__(self,records):
        self.records = [ SimpleSeq(r) for r in list(records) ]
        self.target = self.records[0]
        self.query = self.records[1]
        return

    def lastz(self, mafout=True, program='lastz --chain --strand=plus'):
        # write tempoaray files of records
        with TmpFiles(len(self.records)) as t:
            for i in range(len(self.records)):
                if self.records[i].gapped():
                    self.records[i].removeGaps()
                self.records[i].writefh(t.handles[i], "fasta")
                t.handles[i].close()

            # setup alignment program
            if program.startswith('lastz'):
                cmd = ' '.join([program, t.handles[0].name, t.handles[1].name])
            else:
                sys.exit('Unknown alignment program')

            # run
            child = subprocess.Popen(str(cmd), stdout=subprocess.PIPE, shell=True)

            # parse results (convert)
            blocks = []
            if mafout:
                # use bx-tools to convert to maf
                mafh = t._tmpfile()
                out = bx.align.maf.Writer(mafh)
                for lavBlock in bx.align.lav.Reader(child.stdout):
                    out.write(lavBlock)
                mafh.close()
                # parse MAF file
                mafread = open(mafh.name,'r')
                blocks = []
                for block in MAF.next(mafread):
                    blocks.append(block)
                mafread.close()
            else:
                # parse lav result
                for lavBlock in bx.align.lav.Reader(child.stdout):
                    blocks.append(block)
        return blocks


if __name__ == "__main__":
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC

    a = SeqRecord(Seq("ACGTAGCTAGGCGCGATCTACTCTTACGTAGCTAGGCGCGATCTACTCTTACGGATCTATCTATTACACGGATCTATCTATTAC"), id="one", name="test_one", description="testseq1")
    b = SeqRecord(Seq("AGGGAGACGTTGCGAGACGTAGCTAGGCGCGATCTACTCTTACGGATCTATCTATTACTACCTCTTACGGAT------TATTCC"), id="two", name="test_two", description="testseq2")
    print 'SeqA', SimpleSeq(a)
    print 'SeqB', SimpleSeq(b)
    alnrun = Align([a,b])
    bl = alnrun.lastz(True, 'lastz --chain --strand=plus')
    for block in bl:
        print block
