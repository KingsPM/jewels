#!/usr/bin/env python

import sys

def getField(l,f):
    a = l.split('; ')
    for aa in a:
        if aa.startswith(f):
            return aa.split('"')[1]
    return


if __name__ == "__main__":
    try:
       field = sys.argv[1]
    except:
        print >> sys.stderr, sys.argv[0], '[FIELD] <> IN/OUT'
        sys.exit(1)

    minmax = {}
    counter = 0
    for line in sys.stdin:
        f = line.rstrip().split('\t')
        if len(f) == 9:
            fieldvalue = getField(f[8], field)
            if fieldvalue not in minmax.keys():
                for k in minmax.keys():
                    sys.stdout.write(minmax[k][0] + '\t' + str(minmax[k][1] - 1) + '\t' + str(minmax[k][2]) + '\t' + k + '\t0\t' + minmax[k][3] + '\n')

                counter += 1
                sys.stderr.write(str(counter) + '\r')

                minmax = {}
                minmax[fieldvalue] = [f[0], int(f[3]), int(f[4]), f[6]]
            else:
                try:
                    assert minmax[fieldvalue][3] == f[6]
                except:
                    print >> sys.stderr, "ERROR: not same strand for same feature name (%s)" % (fieldvalue)
                    raise
                if minmax[fieldvalue][1] > int(f[3]):
                    minmax[fieldvalue][1] = int(f[3])
                if minmax[fieldvalue][2] < int(f[4]):
                    minmax[fieldvalue][2] = int(f[4])

    for k in minmax.keys():
        sys.stdout.write(minmax[k][0] + '\t' + str(minmax[k][1] - 1) + '\t' + str(minmax[k][2]) + '\t' + k + '\t0\t' + minmax[k][3] + '\n')

    #scaffold_0	protein_coding	CDS	22976	23049	.	+	0	gene_id "ab.gene.s0.1"; transcript_id "ab.mrna.s0.1.2"; exon_number "5"; protein_id "ab.prot.s0.1.2";
