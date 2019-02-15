#!/usr/bin/env python

import sys
import re

seqidParser = re.compile('(\D+)(\d+)(?:\S+)?\s+(\d+)\s+(\d+)')

data = {}
for line in sys.stdin:
    # create sort key
    m = seqidParser.match(line)
    try:
        sortkey = (m.group(1),int(m.group(2)),int(m.group(3)),int(m.group(4)))
    except:
        sys.stderr.write(line)
        raise
    # store
    try:
        data[sortkey].append(line)
    except:
        data[sortkey] = [line]

# write sorted
for k in sorted(data.keys()):
    sys.stdout.write(''.join(data[k]))
