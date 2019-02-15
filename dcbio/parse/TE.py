#!/usr/bin/env python

import sys

__doc__='''
parser for Shaohua Fan's TE annotation format
'''

#TE_percent_divergence	scaffold_name	TE_name	TE_start	TE_end	TE_orientation	gene/UTR_name	gene/UTR_start	gene/UTR_end	gene/UTR_orientation	TE_gene/UTR relationship
#0.0	LG1	Unknown_R=6174	31120187	31120230	-	on.gene.LG1.753	31137449	31173675	+	locates within the 20kb upstream or downstream of gene
#28.3	LG1	L1-Tx1-1_DR20	31120495	31120806	-	on.gene.LG1.753	31137449	31173675	+	locates within the 20kb upstream or downstream of gene

def read(fh):
    te = []
    for line in fh:
        te.append(TE(line))
    return te

class TE:
    def __init__(self,line):
        f = line.rstrip().split('\t')
        self.divergence = float(f[0])
        self.seqname = f[1]
        self.name  = f[2]
        self.start = int(f[3])
        self.end = int(f[4])
        self.strand = f[5]
        self.geneid = f[6]
        self.genestart = int(f[7])
        self.geneend = int(f[8])
        self.genestrand = f[9]
        self.generelationship = f[10]
        return self

