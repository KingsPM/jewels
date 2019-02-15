#!/usr/bin/env python

'''
extremely simple tree vizualisaion
'''

from Bio import Phylo
import sys

for f in sys.argv[1:]:
    tree = Phylo.read(sys.argv[1], "newick")
    Phylo.draw_ascii(tree)
