#!/usr/bin/env python

'''
similar to grep -f but requires a specified field (default is first)
'''

import sys
from optparse import OptionParser

# read options
parser = OptionParser()
parser.add_option("-v", action="store_true",dest="negate",default=False, help="Negate (like \'grep -v\')")
parser.add_option("-f", dest='termfile',help="search term file [MANDATORY]")
parser.add_option("-c", dest='column',default='1', help="column number [1]")
parser.add_option("-d", dest='delim',default='\t', help="delimier [tab]")
parser.add_option("-m", action="store_true",dest='fullmatch',default=False, help="match entire field [DEFAULT: False]")
(options, args) = parser.parse_args()

print >> sys.stderr, options
print >> sys.stderr, args

# prevent multiple standard inputs
if options.termfile == '-' and (len(args) == 0 or args[0] == '-'):
    parser.error('Cannot use multiple standard inputs (STDIN)')

# configure term handle
if options.termfile is None:   # if filename is not given
    parser.error('The -f option is mandatory (otherwise just use GNU grep...)')
elif options.termfile == '-':
    fh = sys.stdin
else:
    fh = open(options.termfile,'r')

# read terms
terms = set([])
for line in fh:
    terms.add(line.split()[0])
fh.close()

print >> sys.stderr, len(terms), "search terms read"

# configure read handle
if len(args) == 0 or args[0] == '-':
    fh = sys.stdin
else:
    fh = open(args[0],'r')

# parse file
column = int(options.column)-1
foundterms = set([])
for line in fh:
    f = line.rstrip().split(options.delim)
    found = False
    if column >= 0:
        to = column + 1
        searchthis = f[column:to]
    else:
        searchthis = f
    if options.fullmatch:
        for i in searchthis:
            if i in terms:
                found = True
                foundterms.add(i)
                break
    else:
        for t in terms:
            if any([ t in x for x in searchthis]):
                found = True
                foundterms.add(t)
                break
    # XOR
    if found != options.negate:
        sys.stdout.write(line)

if options.fullmatch:
    print >> sys.stderr, "\nMissed terms:\n-------------"
    for t in list(terms.difference(foundterms)):
       print >> sys.stderr, t

