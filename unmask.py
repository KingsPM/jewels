#1/usr/bin/env python

'''removes softmasking from fasta file'''

import sys

for line in sys.stdin:
    if line.startswith('>'):
        sys.stdout.write(line)
    else:
        sys.stdout.write(line.upper())
