#!/usr/bin/env python

"""Matches the QiaSeq linker sequence and
creates a live histogram of UMI wobble

Just pipe in a stream of uncompressed DNA sequences
eg. from a FASTQ file with `zcat R2.fastq.gz | awk '(NR%4==2) { print }'`

a self updating histogram will show
terminates after 2s once stream ends
"""

import sys
import time
from collections import Counter
from curses import wrapper
import regex

__author__ = "David Brawand"
__copyright__ = "KCH"
__license__ = "GPLv3"
__version__ = "0.0.1"

MAX_LEN = 30
MAX_WIDTH = 50.0
LINKER = regex.compile(
    '(.{1,'+str(MAX_LEN)+'})' +
    '(?:[AN][TN]{2}[GN]{2}[AN][GN][TN][CN]{2}[TN]){s<=2,d<=1,i<=1}')
HIST_LINE = '{:>3} |{:<50}| {:8d} ({:.3f})'
HIST_REFRESH = 1000
NOMATCH_KEY = '--'


def main(stdscr):
    """count UMI legths and display histogram"""
    stdscr.clear()

    counter = Counter()
    for i, line in enumerate(sys.stdin):
        matched = LINKER.match(line)
        if matched:
            umi_length = len(matched.group(1))
            counter[umi_length] += 1
        else:
            counter[NOMATCH_KEY] += 1
        # print histogram every 1000 lines
        if i % HIST_REFRESH == 0:
            update_histogram(stdscr, counter)
    update_histogram(stdscr, counter)
    time.sleep(2)


def update_histogram(stdscr, counter):
    """updates the histogram on screen"""
    unit = MAX_WIDTH / float(max(counter.values()))
    total_reads = sum(counter.values())

    for linenumber, i in enumerate([NOMATCH_KEY] + list(range(1, MAX_LEN))):
        stdscr.addstr(linenumber, 0,
                      HIST_LINE.format(i, '*'*int(counter[i] * unit),
                                       counter[i],
                                       counter[i] / float(total_reads)))
    stdscr.refresh()


if __name__ == "__main__":
    wrapper(main)
