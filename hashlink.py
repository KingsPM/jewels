#!/usr/bin/env python
'''
Calculates sequence hashes for 2 multi-fasta files and find common sequences
Returns the corresponding sequence IDs (in order if ambiguous)
If no links can be found assignments are made based on size of sequence (if not ambiguous)
'''

import sys
import md5
from Bio import SeqIO

def seqHash(file):
    counter = 0
    fh = open(file)
    seqh = {}  # seq hash
    seql = {}  # seq length
    seqi = {}  # inversed
    for record in SeqIO.parse(fh, "fasta") :
        counter += 1
        sh = md5.new(str(record.seq).upper()).hexdigest()
        sl = len(record)
        # store seq name
        try:
            seqh[sh].append(record.id)
        except:
            seqh[sh] = [ record.id ]
        try:
            seql[sl].append((record.id, sh))
        except:
            seql[sl] = [ (record.id, sh) ]
        try:
            assert not seqi.has_key(record.id)
        except:
            print 'ERROR: sequence names are not unique (%s)' % (record.id)
            raise
        else:
            seqi[record.id] = (sh,sl)
        # report
        sys.stderr.write('\r' + record.id + ' ' + sh + ' ' + str(sl) + '                          \r')
    fh.close()
    sys.stderr.write('\r' + str(counter) + ' sequences in file ' + file + '                          \n')
    # warn if ambiguity in sequences
    for sh in seqh.keys():
        if len(seqh[sh]) > 1:
            print 'WARNING: Same sequences ' + ','.join(seqh[sh])
    return seqh, seql, seqi

if __name__=="__main__":
    if len(sys.argv) != 3:
        print >> sys.stderr, "USAGE: %s <FASTA> <FASTA>" % (sys.argv[0])
        sys.exit(1)

    a, al, ai = seqHash(sys.argv[1])
    b, bl, bi = seqHash(sys.argv[2])

    # link by hash
    linked = set([])
    for k in a.keys():
        if b.has_key(k):
            linked.add(k)
            for i in range(0,len(a[k])):
                sys.stdout.write(a[k][i] + '\t')
                if i < len(b[k]):
                    sys.stdout.write(b[k][i])
                else:
                    sys.stdout.write('???')
                sys.stdout.write('\tHASH\t' + k + '\n')

    # try to link by length (only unambiguous!)
    unlinked = set(a.keys()) - linked
    for s in unlinked:
        if len(a[s])!=1:
            continue
        name = a[s][0]
        l = ai[name][1]
        if len(al[l])==1:
            try:
                assert len(bl[l]) == 1
            except KeyError, AssertionError:
                pass
            else:
                linked.add(al[l][0][1])
                linked.add(bl[l][0][1])
                sys.stdout.write(al[l][0][0] + '\t' + bl[l][0][0] + '\tLENGTH\t' + str(l) + '\n')

    # print unlinked seqIDs
    for k in a.keys():
        if k not in linked:
            for i in range(0,len(a[k])):
                print a[k][i] + '\t\tUNLINKED\t' + k
    for k in b.keys():
        if k not in linked:
            for i in range(0,len(b[k])):
                print '\t' + b[k][i] + '\tUNLINKED\t' + k
