#!/usr/bin/env python

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import *
import rpy2.robjects as robjects

#import rpy2
#print dir(rpy2.robjects.vectors)
import sys


def lv2py(lv):
    # converts a listvector to python list
    ret = {}
    for key, val in zip(robjects.r.names(lv), lv):
        ret[key] = [i for i in val]  # R Vector to List
    return ret


def readCol(fi, col=1):
    fh = open(fi)
    resu = []
    for line in fh:
        f = line.split()
        if len(f) >= col:
            resu.append(f[col-1])
    fh.close()
    return StrVector(resu)


class GO:
    def __init__(self):
        self.r_base = importr('base')
        self.stats = importr('stats')
        self.topgo = importr('topGO')
        # read functions
        try:
            self.r_base.source('~/biosrc/R/go.R')
            self.goenrichment = robjects.r['testGOquick']
        except:
            print >> sys.stderr, "ERROR: Make sure the go.R file is in an accessible path"
            raise
        # Data
        self.gg = None
        self.fg = None
        self.bg = None

    def gene2go(self, fi):
        # read file and create gene2go list
        fh = open(fi)
        g2g = {}
        for line in fh:
            f = line.split()
            try:
                g = f[1].split(',')
            except:
                g = []
            g2g[f[0]] = StrVector(g)
        fh.close()
        self.gg = ListVector(g2g)

    def testEnrichment(self, ontologies=['MF', 'BP', 'CC']):
        self.resu = self.goenrichment(self.gg, self.bg, self.fg, ontologies)
        return self.resu

    def getSig(self, pval=0.05, adj=False):
        pcol = 8 if adj else 6
        # make list
        resu = lv2py(self.resu)
        ret = {}
        # for each ontology
        for k in resu.keys():
            ret[k] = []
            try:
                rows = max([i for i in range(len(resu[k][pcol])) if float(resu[k][pcol][i]) <= pval])
            except:
                pass
            else:
                # aggregate significant rows
                for r in range(rows+1):
                    ret[k].append([])
                    for c in range(len(resu[k])):
                        ret[k][-1].append(resu[k][c][r])
        return ret

    def testEnrichmentNative(self):
        sys.exit('Native topGO calls not implemented')
        '''# build geneList
        geneList <- robjects.r('factor(as.integer(geneNames%in%selList))')
        names(geneList) <- geneNames



        [ i for i in range(len(listVector2py(self.resu)['MF'][pcol])) if float(listVector2py(g.resu)['MF'][pcol][i]) < 0.05 ]

        # run GO analysis
        results = {}
        for ontology in ('MF','BP','CC'):
            godata <- new("topGOdata", description='session', ontology=onto, allGenes=geneList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO = goList)
            gotest.fisher <- runTest(godata, algorithm='classic', statistic='fisher')
            gotest.ks <- runTest(godata, algorithm='classic', statistic='ks')
            goscores <- score(gotest.fisher)
            gonames <- names(goscores)
            # summarize (top100 enriched)
            num_summarize<-min(c(20, length(gonames)))
            resu[[onto]]<-GenTable(godata, elimFisher=gotest.fisher, elimKS=gotest.ks, orderBy="elimFisher", topNodes=num_summarize)
            # adjust p_values
            resu[[onto]]$elimFisher_adjusted<-p.adjust(resu[[onto]]$elimFisher, method = "BH", n = length(goscores))
            resu[[onto]]$elimKS_adjusted<-p.adjust(resu[[onto]]$elimKS, method = "BH", n = length(goscores))

        self.resu = self.goenrichment(self.gg, self.bg, self.fg)
        return self.resu'''

    def foreground(self, fi):
        try:
            open(fi)
        except:
            self.fg = StrVector(fi)
        else:
            self.fg = readCol(fi)
        return

    def background(self, fi):
        try:
            open(fi)
        except:
            self.bg = StrVector(fi)
        else:
            self.bg = readCol(fi)
        return

    def save(self, file='go.RData'):
        robjects.r.assign('go.bg', self.bg)
        robjects.r.assign('go.fg', self.fg)
        robjects.r.assign('go.gg', self.gg)
        robjects.r.assign('go.resu', self.resu)
        self.r_base.save_image(file)
        return

    def ls(self):
        return


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print >> sys.stderr, sys.argv[0], "foreground", "background", "gene2go"
        sys.exit(1)
    g = GO()
    g.foreground(sys.argv[1])
    g.background(sys.argv[2])
    g.gene2go(sys.argv[3])
    g.testEnrichment(['MF', 'BP'])
    g.save('go.RData')
    resu =  g.getSig(0.05, False)  # <0.05 non-adjusted p-values
    for k in resu.keys():
        print k
        for r in resu[k]:
            rr = tuple(r[:6] + [ float(r[6]),float(r[7]) ] + r[8:])
            #print len(rr)
            print "\t%12s %50s %4d %3d %.2f %4d %1.3f %1.3f %1.2f %1.2f" % rr
            #print r[:6] + [ float(r[6]),float(r[7]) ] + r[8:]
            #print "\t" + '\t'.join([ str(e) for e in r])
    sys.exit("DONE")

