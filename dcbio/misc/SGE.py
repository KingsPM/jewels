#!/usr/bin/env python

import sys
import re
import os
from subprocess import Popen, PIPE
import getpass
from uuid import uuid4

queue = 'medium_jobs.q'

def callQstat(user=getpass.getuser()):
    qs = Popen("qstat", shell=True, stdout=PIPE).stdout
    qstat = ''
    for l in qs:
        if re.match(r"^\d+", l):
            if user in l:
                qstat += re.sub(r' +', r'\t', l)
    return qstat


def batchSGE(commands, name):
    # generate random file name
    randomfile = 'sgejob' + str(uuid4())[24:] + '.sh'
    # write commands to file
    fh = open(randomfile, 'w')
    print >> fh, '#!/bin/bash'
    for c in commands:
        print >> fh, c
    fh.close()
    os.chmod(randomfile, 0755)
    # launch SGE job
    subm = 'qsub -q '+queue+' -b n -V -N ' + name + ' ' + randomfile
    os.system(subm)
    return randomfile

class SGE:
    def __init__(self, run, name, proc, workdir=''):
        self.name = name
        self.proc = str(proc)
        self.queue = queue
        self.err = 'stderr.' + name
        self.out = 'stdout.' + name
        self.command = run
        # GridEngine needs absolute paths for some reason
        wd = ''
        if workdir.startswith('./'):
            wd = os.getcwd() + workdir[1:]
        # build execstring
        self.execstring = 'qsub -V -b y '
        self.execstring += '-wd ' + wd + ' ' if wd else '-cwd '
        self.execstring += '-q ' + self.queue + ' '
        self.execstring += '-N ' + self.name + ' '
        # paml doesn't need that (all run data will be in workdir)
        # self.execstring += '-o ' + wd + self.out + ' '
        # self.execstring += '-e ' + wd + self.err + ' '
        self.execstring += r"'" + self.command + r"'"
        # jobid
        self.jobid = 0  # 0 is no job submitted yet
        return

    def submit(self, debug=False):
        if debug:
            print >> sys.stderr, "## DEBUG ##", self.execstring
        else:
            fh = os.popen(self.execstring)
            for i in fh.readlines():
                self.jobid = int(i.split(' ')[2])
                return self.jobid
        return 0

    def status(self, qstat):
        for i in qstat.split('\n'):
            f = i.split('\t')
            if int(f[0]) == self.jobid:
                return f[4]
        if self.jobid != 0:
            return ''  # done if jobid but not in qstat
        return None  # returns none if job inexistent


