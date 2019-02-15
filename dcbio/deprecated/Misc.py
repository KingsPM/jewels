#!/usr/bin/env python

'''convenience functions'''

import re

## regex and simple match and return function
def rematch(pattern,inp):
    matcher = re.compile(pattern)
    matches = matcher.match(inp)
    if matches:
        yield matches

