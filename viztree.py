#!/usr/bin/env python

import sys
from Bio import Phylo
from cStringIO import StringIO

if __name__=="__main__":
    counter = 0
    sorting = []
    trees = {}
    name = None
    for line in sys.stdin:
        if line.startswith('('):
            counter += 1
            index = name if name else 'Tree' + str(counter)
            trees[index] = Phylo.read(StringIO(line), "newick")
            sorting.append(index)
            name = None
        elif line.startswith('>'):
            name = line[1:].rstrip()
    # plot trees
    headerpre = 5
    headerlen = 80
    print
    for t in sorting:
        print '#'*headerlen
        print '#'*headerpre, t, '#'*(headerlen-headerpre-2-len(t))
        print '#'*headerlen
        print
        Phylo.draw_ascii(trees[t])
        print



