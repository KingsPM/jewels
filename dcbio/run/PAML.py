#!/usr/bin/env python

'''
PAML class
'''

import sys
import os
from Bio.Phylo.PAML import codeml, baseml

class PAML:
    def __init__(self, f, prog):
        # alignment name and bootstrap index (0: full alignment)
        self.name = f[0]
        self.alnum = f[1]
        self.length = f[2]
        self.bootstrap = f[3]
        # alignment tree out_file working_dir seqtype CodonFreq model NSsites
        self.alignment = f[4]
        self.tree = f[5]
        self.out_file = f[6]
        self.working_dir = f[7]
        self.seqtype = f[8]
        self.CodonFreq = f[9]
        self.model = f[10]
        self.NSsites = f[11]

        self.sge = None
        self.prog = prog
        if self.prog == 'codeml':

            self.ml = codeml.Codeml()
            self.ml.alignment = self.working_dir + '/' + self.alignment
            self.ml.tree = self.working_dir + '/' + self.tree
            self.ml.out_file = self.working_dir + '/' + self.out_file
            self.ml.working_dir = self.working_dir

            self.ml.set_options(seqtype=self.seqtype)
            self.ml.set_options(CodonFreq=self.CodonFreq)
            self.ml.set_options(model=self.model)
            self.ml.set_options(NSsites=self.NSsites)
            self.ml.set_options(verbose=1)
            self.ml.set_options(noisy=9)

        elif self.prog == 'baseml':
            self.ml = baseml.Baseml()
            self.ml.alignment = self.alignment
            self.ml.tree = self.tree
            self.ml.out_file = self.out_file
            self.ml.working_dir = self.working_dir
            self.ml.model = 3  # F84 (5:T92 6:TN93)
            self.ml.mgene = 0  # rates
        else:
            print >> sys.stderr, "ERROR: unkown program", self.prog
        return

    def writeCTL(self):
        # write CTL
        try:
            os.mkdir(self.working_dir)
        except OSError:
            pass
        self.ml.ctl_file = self.working_dir + '/' + self.prog + '.ctl'
        self.ml.write_ctl_file()
        return 

    def SGErun(self, debug=False):
        self.writeCTL()
        # return SGE command line
        self.sge = SGE(self.prog + ' ' + self.prog + '.ctl', self.name + '_' + str(self.bootstrap), 1, self.working_dir)
        self.sge.submit(debug)
        # should return job number
        return self.sge

    def SGEcommands(self):
        self.writeCTL()
        # CD to working dir
        commands = ''
        wd = self.working_dir
        if self.working_dir.startswith('./'):
            wd = os.getcwd() + self.working_dir[1:]
        commands += 'cd ' + os.getcwd() + '/' + wd + '\n'
        # run paml
        commands += '/net/isi-cgat/ifs/apps/bio/paml-4.4c/bin/codeml codeml.ctl'
        return commands

