#!/usr/bin/env python

'''convert (flattens) a BED6/GTF file to BED12'''

import sys
#import argparse

def getField(l,f):
    a = l.split('; ')
    for aa in a:
        if aa.startswith(f):
            return aa.split('"')[1]
    return

def writeBedlines(mm, uniq=True, flat=True, bed12=True):
    for k in mm.keys():
        # seqid start stop name score strand
        # tStart tEnd rgb blockCount blockSizes blockStarts
        sys.stdout.write(mm[k][0] + '\t' + str(mm[k][1]) + '\t' + str(mm[k][2]) + '\t' + k + '\t0\t' + mm[k][3]) # BED6
        # BED12 extensions
        if bed12:
            if mm[k][5] is not None and mm[k][6] is not None:
                sys.stdout.write('\t' + str(mm[k][5]))
                sys.stdout.write('\t' + str(mm[k][6]))
            else:
                sys.stdout.write('\t' + str(mm[k][1]))
                sys.stdout.write('\t' + str(mm[k][2]))
            sys.stdout.write('\t' + '128,0,0')

            bl = sorted(mm[k][4])  # sort block tuples
            # unique
            if uniq or flat:
                bl = sorted(list(set(bl)))
            # flattening
            if flat and len(bl) > 1:
                fbl = [list(bl[0])]
                for i in range(1,len(bl)):
                    if fbl[-1][0] + fbl[-1][1] >= bl[i][0]:
                        fbl[-1][1] = max(fbl[-1][0] + fbl[-1][1], bl[i][0] + bl[i][1]) - fbl[-1][0]
                    else:
                        fbl.append(list(bl[i]))
                bl = fbl
            # print blocks
            blocks = zip(*bl)
            sys.stdout.write('\t' + str(len(blocks[0])))
            sys.stdout.write('\t' + ','.join([str(si) for si in list(blocks[1])])) # blockSizes
            sys.stdout.write('\t' + ','.join([str(st-mm[k][1]) for st in list(blocks[0])])) # blockStarts
        sys.stdout.write('\n')
    return

if __name__ == "__main__":

    # bed12
    # unique
    # flatten
    # field
    #### use argparse in future version

    '''modification'''
    '''if from GTF and exon ands CDS are present, thickStart/thickEND get set to the CDS start and end'''

    try:
       field = sys.argv[1]
    except:
        print >> sys.stderr, sys.argv[0], '[FIELD] <> IN/OUT'
        print >> sys.stderr, "NOTE: when inputting a BED file, _name_ should be used, or column number"
        sys.exit(1)


    minmax = {}
    counter = 0
    done = set([])
    for line in sys.stdin:
        cdsStart = None
        cdsEnd = None
        # parse field values
        f = line.rstrip().split('\t')
        if len(f) == 9:
            '''GTF'''
            fieldvalue = getField(f[8], field)
            chrom = f[0]
            blockStart = int(f[3])-1
            blockSize = int(f[4])-blockStart
            strand = f[6]
            blockType = f[2]
            if blockType in ['CDS', 'stop_codon']:
                cdsStart = blockStart
                cdsEnd = blockStart+blockSize
        elif len(f) == 6:
            '''BED6'''
            fieldvalue = f[3] if field == 'name' else f[int(field)-1]
            chrom = f[0]
            blockStart = int(f[1])
            blockSize = int(f[2])-blockStart
            strand = f[5]
            blockType = 'exon'
        elif len(f) == 4:
            '''BED4'''
            fieldvalue = f[3] if field == 'name' else f[int(field)-1]
            chrom = f[0]
            blockStart = int(f[1])
            blockSize = int(f[2])-blockStart
            strand = '.'
            blockType = 'exon'

        # skip start codons as they are included in CDS (stop_codons aren't)
        if blockType in ['start_codon']:
            continue
        # new feature (print old and check if ordered)
        elif fieldvalue not in minmax.keys():
            # sorting check
            try:
                assert fieldvalue not in done
            except:
                sys.stderr.write('ERROR ## input file not sorted by id/field ('+field+')!')
                raise
            else:
                done.add(fieldvalue)

            # write previous
            writeBedlines(minmax)

            # counter
            counter += 1
            sys.stderr.write(str(counter) + '\r')

            # reset and restart
            minmax = {}
            minmax[fieldvalue] = [chrom, blockStart, blockStart+blockSize, strand, [(blockStart,blockSize)], cdsStart, cdsEnd ]
        # extend feature
        else:
            # strand check
            try:
                assert minmax[fieldvalue][3] == strand
            except:
                print >> sys.stderr, "ERROR: not same strand for same feature name (%s)" % (fieldvalue)
                raise

            # adjust ends (thin)
            if minmax[fieldvalue][1] > blockStart:
                minmax[fieldvalue][1] = blockStart
            if minmax[fieldvalue][2] < blockStart+blockSize:
                minmax[fieldvalue][2] = blockStart+blockSize

            # add/adjust ends (thick)
            if blockType in ['CDS', 'stop_codon']:
                if minmax[fieldvalue][5] is None or minmax[fieldvalue][5] > cdsStart:
                    minmax[fieldvalue][5] = cdsStart
                if minmax[fieldvalue][6] is None or minmax[fieldvalue][6] < cdsEnd:
                    minmax[fieldvalue][6] = cdsEnd

            # add block
            minmax[fieldvalue][4].append((blockStart,blockSize))

        #print >> sys.stderr, '###', line,
        #print >> sys.stderr, '###', minmax


    writeBedlines(minmax)

    # scaffold_0	protein_coding	CDS	22976	23049	.	+	0	gene_id "ab.gene.s0.1"; transcript_id "ab.mrna.s0.1.2"; exon_number "5"; protein_id "ab.prot.s0.1.2";
    # seqid start stop name score strand tStart tEnd rgb blockCount blockSizes blockStarts
