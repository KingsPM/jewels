#!/usr/bin/env/python python

import sys
from CIGAR import CIGAR
from dcbio.parse.SequenceFile import SequenceFile, Sequence
from dcbio.misc.SGE import SGE


'''Wrapper class for running and parsing lastZ alignments'''
class LastZ:
    def __init__(self):
        self.run = None
        self.result = None
        return

    def setup(self,target,query,parameters,output):
        self.run = LastZrun(query,target,parameters,output)
        return

    def (self):
        pass


'''LastZ result (CIGAR, start, end)'''
class LastZresult:
    def __init__(self):
        self.cigar = None
        self.alignment = None
        self.targetStart = 
        self.targetEnd = 
        self.queryStart = 
        self.queryEnd = 
        return


'''LastZ factory'''
class LastZrun:
    def __init__(self,query,target,output):
        self.prog   = 'lastz'
        self.query  = query
        self.target = target

    # immediate run
    def SGErun(self, debug=False):
        # return SGE command line
        self.sge = SGE(
            ' '.join([self.prog, self.target, self.query,]),
            '_'.join(self.target.name),
            1,
            self.working_dir
            )
        self.sge.submit(debug)
        # should return job number
        return self.sge

    # returns commands for batch processing
    def SGEcommands(self):
        # CD to working dir
        commands = ''
        wd = self.working_dir
        if self.working_dir.startswith('./'):
            wd = os.getcwd() + self.working_dir[1:]
        commands += 'cd ' + os.getcwd() + '/' + wd + '\n'
        # run command
        commands += ' '.join([self.prog, self.target, self.query,])
        return commands




'''test routine and pipeline protoype'''
if __name__=="__main__":
    try:
        sys.argv[1]
    except:
        raise

    with open(sys.argv[1]) as fh:
        for line in 
