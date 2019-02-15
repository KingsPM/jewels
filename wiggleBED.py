#!/usr/bin/env python

__doc__ = '''
Generates coverage graphs from BED annotated regions
- meta/stacked
- flanks
- reads BED and bigwig
'''

import sys
import re
import os
import struct
from optparse import OptionParser
from bx.bbi.bigwig_file import BigWigFile
import time

import dcbio.parse.BEDfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2014, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


# chromsizes from bigwig
def BigWigChrSize(fh):
    # FROM http://seqanswers.com/forums/showthread.php?t=16888
    csize = {}
    #fh = open(os.path.expanduser(bwname), "rb")
    # read magic number to guess endianness
    magic = fh.read(4)
    if magic == '&\xfc\x8f\x88':
        endianness = '<'
    elif magic == '\x88\x8f\xfc&':
        endianness = '>'
    else:
        raise IOError("The file is not in bigwig format")
    # read the header
    (version, zoomLevels, chromosomeTreeOffset, fullDataOffset, fullIndexOffset, fieldCount, \
        definedFieldCount, autoSqlOffset, totalSummaryOffset, uncompressBufSize, reserved) = \
        struct.unpack(endianness + 'HHQQQHHQQIQ', fh.read(60))
    if version < 3:
        raise IOError("Bigwig files version <3 are not supported")
    # go to the data
    fh.seek(chromosomeTreeOffset)
    # read magic again
    magic = fh.read(4)
    if magic == '\x91\x8c\xcax':
        endianness = '<'
    elif magic == 'x\xca\x8c\x91':
        endianness = '>'
    else:
        raise ValueError("Wrong magic for this bigwig data file")
    (blockSize, keySize, valSize, itemCount, reserved) = struct.unpack(endianness + 'IIIQQ', fh.read(28))
    (isLeaf, reserved, count) = struct.unpack(endianness + 'BBH', fh.read(4))
    for n in range(count):
        (key, chromId, chromSize) = struct.unpack(endianness + str(keySize) + 'sII', fh.read(keySize + 2 * 4))
        # we have chrom and size
        csize[key.replace('\x00', '')] = chromSize
    # rewind filehandle
    wigfh.seek(0)
    return csize

# fast running mean (using convolution)
def _runningMeanFast(x, N):
    # get running mean
    return np.convolve(x, np.ones((N,))/N)[(N-1):]


def populate(wigfh,bedfh,window=1,prefix='coverage'):
    chromsizes = BigWigChrSize(wigfh)
    bw = BigWigFile(wigfh)
    annotation = dcbio.parse.BEDfile.BED(bedfh)
    for b in annotation:
        # extreact coverage
        start = b.chromStart - options.flank if options.flank <= b.chromStart else 0
        end = b.chromEnd + options.flank if b.chromEnd + options.flank <= chromsizes[b.chrom] else chromsizes[b.chrom]

        #datapoints = len(b) / int(options.res)
        #b.meta['coverage'] = bw.query(b.chrom, start, end, datapoints)
        b.meta[prefix] = [ x[2] for x in bw.get(b.chrom, start, end) ]
        b.meta['_'.join([prefix,'mean'])] = np.mean(b.meta['coverage'])
        b.meta['_'.join([prefix,'window'])] = window
        if window != 1:
            b.meta[prefix] = _runningMeanFast(b.meta[prefix],window)  # right facing windows
        print >> sys.stderr, "POP", b.meta
    # reset filehandles
    wigfh.seek(0)
    return annotation

## MAIN
if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-w", "--bigwig",     dest="wigfile", metavar="STRING", help="BigWig file")
    parser.add_option("-b", "--bed",        dest="bedfile", metavar="STRING", help="BED file")
    parser.add_option("-f", "--flanking",   dest="flank",   metavar="INT", default=0, help="flanking sequence [0]")
    parser.add_option("-r", "--resolution", dest="res",     metavar="INT", default=1, help="resolution [1]")
    (options, args) = parser.parse_args()

    # open bigwig, get chrom sizes
    try:
        wigfh = open( options.wigfile )
    except:
        if options.wigfile == '-':
            wigfh = sys.stdin
        else:
            raise Exception("WIG file cannot be openend/found")

    # open bedfile or select stdin
    try:
        bedfh = open(options.bedfile)
    except:
        if options.bedfile == '-':
            bedfh = sys.stdin
        else:
            raise Exception("BED file cannot be openend/found")

    # get coverage per BEDline
    coverage = populate(wigfh,bedfh,int(options.res))

    # get coding exons from bed file and intersect


    # assess coverage over coding exons


'''
        plt.plot([ x[0] for x in b.meta['coverage']], [ x[2] for x in b.meta['coverage']], 'r-')
        plt.ylabel('read depth')
        plt.xlabel('coordinate')
        plt.show()
        sys.exit('DEBUG')
'''


'''
read shaohuas data
compare with etienne
make this thing more interactive
'''


