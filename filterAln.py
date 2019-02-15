#!/usr/bin/env python

__doc__ = '''
simple sliding window filtering based on mismatch threshold
removes non-aligned sequence prior to filtering
'''

__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2013, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]

import sys
import dcbio.parse.SimpleAlignment
from shutil import copyfile
from optparse import OptionParser

if __name__=='__main__':
    # USAGE
    usage = ' '.join(["USAGE:",sys.argv[0],"[options]","ALIGNMENT_FILE"])
    parser = OptionParser(usage=usage)
    parser.add_option("-t", dest="threshold", metavar="INT", default=6,  help="maximum number of mismatched bases in window [6]")
    parser.add_option("-w", dest="window",    metavar="INT", default=12, help="sliding window size [12]")
    parser.add_option("-s", dest="slide",     metavar="INT", default=3,  help="sliding distance [3]")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print >> sys.stderr, usage
        sys.exit(1)

    # read alignment and rename file
    aln = SimpleAlignment.Aln()
    suffix = args[0].split('.')[-1]
    with open(args[0]) as fh:
        if 'fa' in suffix:
            aln.fasta(fh)
            outputformat = 'fa'
        elif 'phy' in suffix:
            aln.phylip(fh)
            outputformat = 'phy'
        elif 'aln' in suffix:
            aln.clustalw(fh)
            outputformat = 'aln'
        else:
            print >> sys.stderr, "(!) unkown file extension (fa,phy,aln)"
            sys.exit(1)

    # remove not aligned sequence (gaps), filter alignment and display
    aln.removeGappedCols()
    filtered, removed = aln.filter(int(options.threshold), int(options.window), int(options.slide))
    aln.display(length=90, color=True, codons=True, blink=removed)

    # write filtered alignment and report result
    outfile = args[0]+'.filtered'
    fh = open(outfile,'w')
    filtered.writefh(fh, outputformat)
    fh.close()
    print outfile, aln.inlength, filtered.inlength, "=" if filtered.inlength == aln.inlength else '-'
