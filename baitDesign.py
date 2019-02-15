#!/usr/bin/env python

__doc__ = '''
picks regions for Bait capture according to specifications (tab-delimited file)
--> force excludes introns for secondary annotations
'''
__author__ = "David Brawand, PhD - Viapath | King's College Hospital"
__copyright__ = "Copyright 2014, David Brawand"
__credits__ = []
__license__ = "MIT"
__version__ = "0.5"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Development"  # ["Prototype", "Development",  "Production"]


import sys
from optparse import OptionParser
from dcbio.parse.BEDfile import BED

class Spec:
    def __init__(self,line):
        f = line.rstrip().split('\t')
        f = [ x.strip('"') for x in f ]
        self.line = line.rstrip()
        self.coordinates = []
        self.region = None

        # sanity check
        try:
           self.number = int(f[0])
        except:
            raise Exception("query not numbered")
        self.geneName = f[1].rstrip() if f[1].startswith('rs') else f[1].rstrip().upper()
        self.category = f[2]
        self.utr = False
        if f[3].upper().startswith('Y'):
            self.utr = True

        # assume if nothing then inlcude introns
        self.introns = False
        try:
            self.flanks = int(f[4])
        except ValueError:
            self.flanks = 30
            self.introns = True
        except IndexError:
            print '???', line.split()
            raise

        # upstream regions
        try:
            self.upstream = abs(int(f[5]))
        except ValueError:
            self.upstream = 0
        except IndexError:
            self.upstream = 0

        # gap closing
        try:
            self.gapclosing = abs(int(f[6]))
        except ValueError:
            self.gapclosing = 0
        except IndexError:
            self.gapclosing = 0

        # secondary annotations (introns are not included)
        try:
            self.secondary = map(int,f[7].split(';'))
        except ValueError:
            self.secondary = []
        except IndexError:
            self.secondary = []


        # region
        try:
            self.region = f[8]
        except IndexError:
            pass
        else:
            if f[8].startswith('chr'):
                cc = [ x.split('-') for x in f[8].replace(',','').split(':')]
                self.coordinates.append([ cc[0][0], int(cc[1][0])-self.flanks, int(cc[1][1])+self.flanks, self.geneName ])

        return

    def __str__(self):
        return self.line

    def __repr__(self):
        reprData = [ self.geneName,
                    "UTR="+str(self.utr),
                    "INTRONS="+str(self.introns),
                    "FLANKS="+str(self.flanks),
                    "PROMOTER="+str(self.upstream),
                    "CLOSEGAPS="+str(self.gapclosing),
                    "SECONDARY="+';'.join(map(str,self.secondary)) ]
        return '\t'.join(reprData)

    def flattenAndCloseGaps(self):
        self.coordinates = sorted(self.coordinates)
        flattened = [ self.coordinates[0] ]
        for c in self.coordinates[1:]:
            if c[1] <= flattened[-1][2] + self.gapclosing:
                # check sorting and uniformity
                try:
                    assert flattened[-1][0] == c[0]  # same chromosome
                    assert flattened[-1][1] <= c[1]  # same or greater start
                    assert flattened[-1][3] == c[3]  # same name
                except:
                    raise
                # extend
                flattened[-1][2] = max(flattened[-1][2],c[2])
            else:
                # new segment
                flattened.append(c)
        # apply
        self.coordinates = flattened
        return

    def BED(self):
        lines = []
        for c in self.coordinates:
            lines.append("\t".join(map(str,c)))
        return "\n".join(lines)

def readSpecs(specfi):
    regions = []
    lastnumber = 0
    with open(specfi) as fh:
        for line in fh:
            try:
                entrynumber = int(line.split()[0])
            except IndexError:
                #print >> sys.stderr, "##### EMPTY #", line.rstrip()
                continue
            except ValueError:
                #print >> sys.stderr, "## SKIPPING #", line.rstrip()
                continue
            else:
                #print >> sys.stderr, "#### TARGET #", line.rstrip()
                pass

            # test numbering (just for sanity)
            try:
                assert entrynumber - lastnumber == 1
            except AssertionError:
                raise Exception("Queries not numbered continuously")
            else:
                lastnumber = entrynumber
                regions.append(Spec(line))
                print >> sys.stderr, repr(regions[-1])
    return regions

def readSynonyms(fi):
    synonyms = {}
    ambiguous = set([])
    with open(fi) as fh:
        for line in fh:
            f = line.split()
            if f[1].startswith('-'):
                continue
            else:
                for s in f[1].split('|'):
                    try:
                        synonyms[s]
                    except KeyError:
                        synonyms[s] = f[0]
                    else:
                        ambiguous.add(s)
    # remove ambiguous synonyms
    for a in ambiguous:
        del synonyms[a]
    return synonyms

if __name__=="__main__":

    usage = "usage: %prog [options] <ANNOTATION1> [ANNOTATION2] ..."
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--specs", dest="specs", metavar="FILE", help="Specification file (tab from excel)")
    parser.add_option("-y", "--synonyms", dest="synonyms", metavar="FILE", help="synonyms (from UCSC gene_info")
    (options, args) = parser.parse_args()

    # read reference BED12 files
    refs = []
    refindex = []
    for refFile in args:
        print >> sys.stderr, "READING annotations from %s" % refFile
        with open(refFile) as fh:
            refs.append(BED(fh))
            refindex.append(refs[-1].getDictList())
        print >> sys.stderr, "READ gene annotations from %s" % refFile

    # read synonyms
    synonyms = {}
    if options.synonyms:
        print >> sys.stderr, "READING synonmyms"
        synonyms = readSynonyms(options.synonyms)
        print >> sys.stderr, "READ gene symbol synonyms"

    # read specification
    # number, gene, category, utr3, utr5, exonflanks, upstream, description
    specs = readSpecs(options.specs)

    # pick regions and output BED (covered regions)
    for s in specs:
        #print >> sys.stderr, "##### FETCH ##", s, '\r',
        # skip the big regions
        if s.region:
            assert s.coordinates
        else:
            # get gene info (primary and secondary)
            genes = None
            secondary = []

            # check if gene symbol exists as such or as synonym
            if s.geneName in refindex[0].keys():
                queryname = s.geneName
            else:
                try:
                    queryname = synonyms[s.geneName]
                except KeyError:
                    print >> sys.stderr, "\n### NOT FOUND ###", str(s.geneName)
                    raise
                else:
                    s.geneName += '|' + queryname

            # select primary
            genes = refindex[0][queryname]
            # select secondary (no introns added)
            for r in [ refindex[i] for i in s.secondary ]:
                try:
                    secondary += r[queryname]
                except:
                    pass

            # fetch primary gene
            try:
                assert genes
            except AssertionError:
                print >> sys.stderr, "\n### NOT FOUND ###", str(s)
                #raise Exception("gene not found")
            else:
                for gene in genes:
                    if s.introns:
                        # full gene (plus flanks plus upstream)
                        if s.utr:
                            if gene.strand == "+":
                                s.coordinates.append([ gene.chrom, gene.chromStart-s.upstream, gene.chromEnd, s.geneName ])
                            elif gene.strand == "-":
                                s.coordinates.append([ gene.chrom, gene.chromStart, gene.chromEnd+s.upstream, s.geneName ])
                            else:
                                raise Exception('Strand error')
                        else:
                            if gene.strand == "+":
                                s.coordinates.append([ gene.chrom, gene.thickStart-s.upstream, gene.thickEnd, s.geneName ])
                            elif gene.strand == "-":
                                s.coordinates.append([ gene.chrom, gene.thickStart, gene.thickEnd+s.upstream, s.geneName ])
                            else:
                                raise Exception('Strand error')
                    else:
                        # exons (flanks and upstream)
                        # extract segments and write
                        cds = not s.utr
                        blocks = gene.getBlocks(thick=cds)

                        # extend into introns and upstreams
                        if gene.strand == '+':
                            s.coordinates.append([ gene.chrom, blocks[0][0]-s.upstream, blocks[0][1]+s.flanks, s.geneName ])
                            for b in blocks[1:]:
                                s.coordinates.append([ gene.chrom, b[0]-s.flanks, b[1]+s.flanks, s.geneName ])
                        elif gene.strand == '-':
                            for b in blocks[:-1]:
                                s.coordinates.append([ gene.chrom, b[0]-s.flanks, b[1]+s.flanks, s.geneName ])
                            s.coordinates.append([ gene.chrom, blocks[-1][0]-s.flanks, blocks[-1][1]+s.upstream, s.geneName ])
                        else:
                            raise Exception('Strand error')

            # fetch secondary genes (don't add any introns)
            try:
                assert secondary
            except AssertionError:
                pass  # no secondary annotation found
            else:
                for gene in secondary:
                    # exons (flanks and upstream)
                    # extract segments and write
                    cds = not s.utr
                    blocks = gene.getBlocks(thick=cds)

                    # extend into introns and upstreams
                    if gene.strand == '+':
                        s.coordinates.append([ gene.chrom, blocks[0][0]-s.upstream, blocks[0][1]+s.flanks, s.geneName ])
                        for b in blocks[1:]:
                            s.coordinates.append([ gene.chrom, b[0]-s.flanks, b[1]+s.flanks, s.geneName ])
                    elif gene.strand == '-':
                        for b in blocks[:-1]:
                            s.coordinates.append([ gene.chrom, b[0]-s.flanks, b[1]+s.flanks, s.geneName ])
                        s.coordinates.append([ gene.chrom, blocks[-1][0]-s.flanks, blocks[-1][1]+s.upstream, s.geneName ])
                    else:
                        raise Exception('Strand error')

        # print blocks
        try:
            s.flattenAndCloseGaps()  # flatten coordinates for more concise output
            print s.BED()
        except:
            print >> sys.stderr, s.coordinates
            raise
