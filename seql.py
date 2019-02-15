#1/usr/bin/env python

import sys

current = None
for line in sys.stdin:
    if line.startswith('>'):
        if current:
            print count, current
        # new query
        current = line[1:-1]
        count = 0
    else:
        count += len(line.rstrip())
# printing
print count, current
