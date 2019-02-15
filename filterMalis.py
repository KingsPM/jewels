#!/usr/bin/env python


'''simple window filtering an alignment in the terminal'''

import sys
import dcbio.parse.SimpleAlignment
from shutil import copyfile

threshold = 2
window = 4
stepping = 1

if len(sys.argv) > 1:
    if sys.argv[1].startswith('c'):
        # codon
        threshold = 5
        window = 12
        stepping = 3
    else:
        sys.exit('Use c for codon filtering (or none for defaults)')

if __name__=='__main__':
    current = []
    for line in sys.stdin:
        if len(current) > 0:
            if current[-1].split()[0] == line.split()[0]:
                current.append(line)
            else:
                aln = SimpleAlignment.Aln()
                aln.malis(current)
                #aln.display()
                aln.removeGappedCols()
                if aln.inlength > 0:
                    filtered, removed = aln.filter(threshold, window, stepping)
                    aln.display(length=90, color=True, codons=True, blink=removed)
                    if filtered.inlength > 0:
                        filtered.writefh(sys.stdout,'malis')
                current = [ line ]
        else:
            current.append(line)

    if len(current) > 1:
        aln = SimpleAlignment.Aln()
        aln.malis(current)
        #aln.display()
        aln.removeGappedCols()
        if aln.inlength > 0:
            filtered, removed = aln.filter(threshold, window, stepping)
            #aln.display(length=90, color=True, codons=True, blink=removed)
            if filtered.inlength > 0:
                filtered.writefh(sys.stdout,'malis')
        current = []

