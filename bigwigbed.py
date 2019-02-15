#!/usr/bin/env python



__doc__ = """
Returns the mean/average coverage in a bigWig file over bed intervals
USAGE: %prog bigwig_file.bw < bed_file.bed
"""

from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile

import numpy as np
import sys

counter = 0
bw = BigWigFile( open( sys.argv[1] ) )
for interval in GenomicIntervalReader( sys.stdin ):
    counter += 1
    cov = [ x[2] for x in bw.get(interval.chrom, interval.start, interval.end)]
    if counter % 100 == 0:
        print >> sys.stderr, counter, '\r',
    print "\t".join(interval.fields + [str(np.mean(cov))])
