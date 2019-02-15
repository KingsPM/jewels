#!/usr/bin/env python

'''live word count'''
'''standard option is counting lines (word count no implemented yet)'''

import sys
import re
from optparse import OptionParser

# parse options
parser = OptionParser()
parser.add_option("-l", "--lines", action="store_true", dest="lines", default=True, help="count lines")
parser.add_option("-w", "--words", action="store_true", dest="words", default=False, help="count words")
parser.add_option("-s","--stepping", dest="step", default=100, metavar="INT")
(options, args) = parser.parse_args()

# open filehandle
if len(args)==0:
    fh = sys.stdin
else:
    fh = open(args[0],'r')

# count
lines, words = 0, 0
for line in fh:
    if options.lines:
        lines += 1
        if lines % options.step == 0:
            sys.stderr.write('\r'+str(lines)+' Lines')
    if options.words:
        words += len(re.split(r'[^0-9A-Za-z]+',line))
        if words % options.step == 0:
            sys.stderr.write('\r'+str(words)+' Words')
if options.lines:
    sys.stderr.write('\r'+str(lines)+' Lines\n')
if options.words:
    sys.stderr.write('\r'+str(words)+' Words\n')

# close filehandle
try:
    fh.close()
except:
    pass
