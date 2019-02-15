#!/usr/bin/env python

__doc__ = '''
select specific blocks from a MAF alignment
'''

import sys
import re
import dcbio.parse.MAF
from optparse import OptionParser

__author__ = "David Brawand, PhD - University of Oxford"
__copyright__ = "Copyright 2014, David Brawand"
__credits__ = []
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "David Brawand"
__email__ = "david.brawand@dpag.ox.ac.uk"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


def parseSpeciesSeqID(s, sep):
    items = set(s.split(','))
    r1, r2 = [], []
    for i in items:
        if sep in i:
            r2.append(i)
        else:
            r1.append(i)
    return r1, r2


if __name__=="__main__":
    # read requirements/includes
    usage = "usage: %prog [options] <MAF>"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", dest="require", metavar="STRING", help="REQUIRE comma separated species(:chromosome) list")
    parser.add_option("-i", dest="include", metavar="STRING", help="INCLUDE comma separated species(:chromosome) list")
    parser.add_option("-s", dest="sep", metavar="CHAR", default=':', help="separator for species and chromsome [:]")
    (options, args) = parser.parse_args()
    if options.require:
        require, require2 = parseSpeciesSeqID(options.require, options.sep)
    else:
        require, require2 = [], []
    if options.include:
        include, include2 = parseSpeciesSeqID(options.include, options.sep)
    elif not options.require:
        parser.print_help()
        sys.exit(1)
    else:
        include, include2 = [], []

    # select read handle
    try:
        readhandle = open(args[0],'r')
    except:
        if len(args)==0 or args[0]=='-':
            sys.stderr.write('Using STDIN for MAF alignment\n')
            readhandle = sys.stdin
        else:
            sys.stderr.write('ERROR: Cannot open file ('+args[0]+')')
            raise
    # select write handle and print header
    writehandle = sys.stdout
    print >> writehandle, '##maf version=1\n#\n'

    counter, success = 0, 0
    # read file
    for block in MAF.next(readhandle,':'):
        # counter
        counter += 1
        # test if requirements for species and sequence are met
        if require and not set(block.species()).issuperset(require):
            # test species requirement failed
            continue
        if require2 and not set(block.speciesSequence(options.sep)).issuperset(require2):
            # test sequence requirement
            continue
        # strip from non-requested stuff
        include_species = block.IDs(True,False) if not include else include
        include_chromosomes = block.IDs(True,True) if not include2 else include2
        block.keep_sequences(spec=include_species, chrom=include_chromosomes)

        if len(block.seqs) > 1:
            success += 1
            block.remove_empty_columns()
            print repr(block)
        else:
            continue
    print >> sys.stderr, "DONE"

    sys.exit(0)
