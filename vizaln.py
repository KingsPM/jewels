#!/usr/bin/env python

__doc__ = '''
visualizes an alignment file in the terminal in color
'''

__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2013, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Production"  # ["Prototype", "Development",  "Production"]

import sys
import dcbio.parse.SimpleAlignment

if __name__=='__main__':
    aln = SimpleAlignment.Aln()

    if len(sys.argv)<2:
        print >> sys.stderr, "USAGE (FILE):", sys.argv[0], "[ALIGNMENT_FILE]"
        print >> sys.stderr, "USAGE (PIPE):", sys.argv[0], "[FORMAT] < ALIGNMENT_FILE"
        sys.exit(1)

    # open file/handle
    try:
        fh = open(sys.argv[1])
    except:
        fh = sys.stdin
        suffix = sys.argv[1]
    else:
        suffix = sys.argv[1].split('.')[-1]

    # parse
    if 'fa' in suffix:
        aln.fasta(fh)
    elif 'aln' in suffix:
        aln.clustalw(fh)
    elif 'phy' in suffix:
        aln.phylip(fh)
    elif 'maf' in suffix:
        aln.maf(fh)
    else:
        print >> sys.stderr, "(!) UNKOWN FILE EXTENSION:", suffix
        sys.exit(1)
    aln.display()  #display(self, length=90, color=True, codons=True, blink=None, fh=sys.stderr)
