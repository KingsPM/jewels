#!/usr/bin/env python

import sys
from tabulate import tabulate

data = []

for line in sys.stdin:
    data.append(line.split())

print tabulate(data)