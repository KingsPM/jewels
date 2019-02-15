#!/usr/bin/env python

#LG1     alignAssembly-pasa_niloticus    cDNA_match      903473  903908  .       -       .       ID=chain_271;Target=asmbl_271 1 436 +
#LG1     alignAssembly-pasa_niloticus    cDNA_match      906310  906957  .       +       .       ID=chain_273;Target=asmbl_273 183 830 +

import sys
import os
import copy
import re

from numpy import std, median, mean

##UNIVERSAL PARSER##
def detect(infile):
    # autodetect filetype
    ext = None
    testFile = open(infile,'r')
    for line in testFile:
        if line.startswith('#'):
            continue
        else:
            f = line.split('\t')
            if "=" in f[8] and not ' "' in f[8]:
                ext = 'gff'
                break
            elif ' "' in f[8] and not '=' in f[8]:
                ext = 'gtf'
                break
    testFile.close()
    if not ext:
        printWarning('Format unknown, assuming GFF3 formatting for file %s' % infile)
        ext = 'gff'
    sys.stderr.write('Format is ' + ext + '\n')
    return ext

def parse(infile,ext='gff'):
    sys.stderr.write('Reading ' + infile + '\n')
    exons = []
    # read by line
    linecount = 0
    fh = open(infile,'r')
    for line in fh:
        linecount += 1
        if linecount % 1000 == 0:
            sys.stderr.write("\r" + str(linecount))
        if line.startswith('#') or len(line.rstrip()) == 0:
            continue
        if ext == 'gtf':
            exons.append(GTF(line))
        elif ext == 'gff':
            exons.append(Generic(line))
        else:
            sys.exit('unknown format')
    fh.close()
    sys.stderr.write('\rFound ' + str(len(exons)) + ' annotation objects\n' )
    return exons

def ident(exons):
    unident = 0 
    for e in exons:
        if not e.hasAttribute('ID'):
            unident += 1
    # set unique IDs if there are none
    if unident > 0:
        for e in exons:
            if not e.hasAttribute('ID'):
                e.setUniqueIdent("U")
    return exons, unident

def link(exons):
    orphans = 0
    for e in exons:
        if not e.hasAttribute('Parent'):
            orphans += 1
    # link method (GFF)
    if orphans != len(exons):
        # make IDs fully unique (CDS)
        uniqID = {}
        parents = set([])
        redun = set([])
        for e in exons:
            # set an id if unavailable
            if not e.hasAttribute('ID'):
                e.setUniqueIdent('UNIQ')
            # determine uniqe ids
            if not uniqID.has_key(e.getAttribute('ID')):
                uniqID[e.getAttribute('ID')] = 0
            if e._strand == '-':
                uniqID[e.getAttribute('ID')] -= 1
            else:
                uniqID[e.getAttribute('ID')] += 1
            # add to redundantID
            if abs(uniqID[e.getAttribute('ID')]) > 1:
                redun.add(e.getAttribute('ID'))
            # make parent set
            if e.hasAttribute('Parent'):
                parents.add(e.getAttribute('Parent'))

        #Set negatives to zero -> absolute values number the exons 
        for k in uniqID.keys():
            if uniqID[k] < 0:
                uniqID[k] = -1

        for e in reversed(exons):
            if e.hasAttribute('ID') and e.getAttribute('ID') in redun:
                if e.hasAttribute('ID') in parents:
                    suicide('parent with non unique ID found (' + e.getAttribute('ID') + ')')
                else:
                    # renumber it
                    newID = e.getAttribute('ID') + '.' + e._type.lower() + str(abs(uniqID[e.getAttribute('ID')]))
                    uniqID[e.getAttribute('ID')] -= 1
                    e.setAttribute('ID='+newID,warn=False)

        # determine if PASA format with gene clusters
        isPasa = False
        isAugustus = False
        if 'AUGUSTUS' in exons[0]._source:
            isAugustus = True
        elif 'pasa' in exons[0]._source:
            isPasa = True
        # PASA-like
        if isPasa:
            sys.stderr.write('Linking (pasa GFF3)')
            # merge gene clusters and modify childs
            geneless = []
            # index mRNA
            mRNAindex = {}
            for e in exons:
                # exclude genes as they will be rebuilt
                if not e._type == 'gene':
                    if e._type == 'mRNA':
                        # truncate parentID
                        parentID = e.getAttribute('Parent')
                        gene_cluster = parentID[:parentID.rfind('.')]
                        e.setAttribute('Parent='+gene_cluster)
                        if not mRNAindex.has_key(gene_cluster):
                            mRNAindex[gene_cluster] = []
                        mRNAindex[gene_cluster].append(e)
                    geneless.append(e)
            # build genes and add to "geneless" set
            for t in mRNAindex.keys():
               geneless.append(Parent(mRNAindex[t],'gene'))
               geneless[-1].setAttribute('ID='+t)
            # sort and replace
            exons = sorted(geneless)
        # augustus (transcript -> mRNA)
        elif isAugustus:
            sys.stderr.write('Linking (augustus GFF3)')
            for e in exons:
                if e._type == 'transcript':
                    e._type = 'mRNA'
        # generic GFF3
        else:
            sys.stderr.write('Linking (generic GFF3)')
        exons = linkAnnotation(exons)
    else:
        # get attributes
        attribs = set()
        for e in exons:
            for a in e._attributes.keys():
                attribs.add(a)
        
        # ensembl (does not have gene and transcripts => add)
        if 'gene_id' in attribs and 'transcript_id' in attribs:
            sys.stderr.write('Linking (ENSEMBL)')
            # index and collect exons per gene and transcript (handles split genes)
            indexed = {}
            for e in exons:
                g = e.getAttribute('gene_id') 
                t = e.getAttribute('transcript_id')
                h = str(e._seqid) + '|' + e.getStrand() # strand seqid hash (used for splitting superstructures which may happen during the liftOver)
                if not indexed.has_key(g):
                    indexed[g] = {}
                if not indexed[g].has_key(h):
                    indexed[g][h] = {}
                if not indexed[g][h].has_key(t):
                    indexed[g][h][t] = []
                indexed[g][h][t].append(e)
            #build transcripts and genes (and link childs)
            added = []
            for g in indexed.keys():
                # split
                subnumbering = 1
                for h in indexed[g].keys():
                    builtTranscripts = []
                    for t in indexed[g][h].keys():
                        # update subID if exons
                        for child in indexed[g][h][t]:
                            if len(indexed[g].keys()) > 1:
                                child.addSubId(str(subnumbering),['transcript_id','gene_id'])
                            child.setAttribute('ID=' + child.getAttribute('transcript_id') + '.' + child._type + child.getAttribute('exon_number'))
                        # build transcript and set ID
                        builtTranscripts.append(Parent(indexed[g][h][t],'mRNA'))
                        builtTranscripts[-1].setAttribute('ID=' + builtTranscripts[-1].getAttribute('transcript_id'))
                        # link
                        for child in indexed[g][h][t]:
                            builtTranscripts[-1].addChild(child)
                    added.extend(builtTranscripts)
                    # build gene
                    added.append(Parent(builtTranscripts,'gene'))
                    added[-1].setAttribute('ID=' + builtTranscripts[-1].getAttribute('gene_id'))
                    for child in builtTranscripts:
                        added[-1].addChild(child)
                    subnumbering += 1
            exons.extend(added)
        # evigan-like (single mRNA per gene)
        elif 'GenePrediction' in attribs:
            sys.stderr.write('Linking (EVIGAN)\n')
            # set ID
            for e in exons:
                if not e.hasAttribute('ID',True):
                    e.setUniqueIdent("U")
            # create parent index
            index = {}
            for e in exons:
                if e._type in set(['mRNA','gene']):
                    name = e.getAttribute('GenePrediction')
                    if not index.has_key(e._type):
                        index[e._type] = {}
                    if not index[e._type].has_key(name):
                        index[e._type][name] = []
                    index[e._type][name].append(e)
            # do the linking
            for e in exons:
                name = e.getAttribute('GenePrediction')
                if e._type == 'gene':
                    pass
                elif e._type == 'mRNA':
                    index['gene'][name][0].addChild(e,False,True) # no check, add parent
                else:
                    index['mRNA'][name][0].addChild(e,False,True) # no check, add parent
                    
        else:
            printWarning('cannot build hierarchy (unknown structure)')
        exons = linkAnnotation(exons,True) # do a linking rebuild to avoid duplicates (this is actually crappy, should have used a set for childs) -> returns toplevel
    return exons

def writefh(exons,fh=sys.stdout,maxdepth=9,fileFormat='GFF'):
    for e in sorted(exons):
        e.write(fh,True,True,maxdepth,fileFormat)
    return

##GLOBAL##
def linkAnnotation(l,rebuild=False):
    # make toplevel index
    index = {}
    for e in l:
        if not e.hasAttribute('ID'):
            print >> sys.stderr, "## WARNING ## no ID for", repr(e)
        ident = e.getAttribute('ID')
        try:
            assert not index.has_key(ident)
        except:
            print '#1',repr(index[ident])
            print '#2',repr(e)
            raise
            suicide('non unique identifyers in file (' + e.getAttribute('ID') + ')')
        index[ident] = e
    # cleanup all child associations
    if rebuild:
        for e in l:
            e._childs = []
    # add childs (no check needed)
    toplevel = []
    for e in l:
        if e.hasAttribute('Parent'):
            parentID = e.getAttribute('Parent')
            index[parentID].addChild(e,False)
        else:
            toplevel.append(e)
    # return toplevel
    return toplevel

def genesFromExons(l):
    # make toplevel index
    index = {}
    for e in l:
        ident = e.getAttribute('ID')
        try:
            assert not index.has_key(ident)
        except:
            suicide('non unique identifyers in file')
        index[ident] = e
    # add childs (no check needed)
    for e in l:
        if e.hasAttribute('Parent'):
            parentID = e.getAttribute('Parent')
            index[parentID].addChild(e,False)
        else:
            toplevel.append(e)
    # return toplevel
    return toplevel

def countAnnotation(l):
    lset = set(l)
    lcount = {}
    levelIndex = {}
    for a in l:
        # determine level
        level = 0
        b = a
        while b._parent is not None:
            level += 1
            b = b._parent
        # count
        if not lcount.has_key(level):
            lcount[level] = 0
        lcount[level] += 1
        levelIndex[a.getAttribute('ID')] = level

    lchilds = { 'data':{}, 'median':{}, 'mean':{}, 'std':{} } # median, mean, stdev
    lvls = set()
    for a in l:
        level = levelIndex[a.getAttribute('ID')]
        lvls.add(level)
        if level < len(lcount) - 1:
            # calc stat
            level = str(level)
            stat = float(len(set(a._childs) & lset))
            try:
                lchilds['data'][level].append(stat)
            except:
                lchilds['data'][level] = [ stat ]
    for l in sorted(lchilds['data'].keys()):
        lchilds['median'][str(l)] = median(lchilds['data'][str(l)])
        lchilds['mean'][str(l)] = mean(lchilds['data'][str(l)])
        lchilds['std'][str(l)] = std(lchilds['data'][str(l)])
    return lcount, lchilds

def keepTypes(l,t,v=1):
    keep = set(t)
    kept = []
    for e in l:
        if e._type in keep:
            kept.append(e)
    if v > 0:
        print >> sys.stderr, "\n# KEPT", len(kept), "OUT OF", len(l), "ANNOTATIONS"
    return kept

def expandHierarchy(p):
    dh = set(p)
    while len(p) > 0:
        siblings = []
        for g in p:
            siblings.extend(g._childs) 
        dh.update(set(siblings))
        p = siblings
    return list(dh)

def covered(features):
    sf = sorted(features)
    cov = [ [sf[0]._start, sf[0]._end] ]
    for i in range(1,len(sf)):
        if sf[i]._start <= cov[-1][1] + 1:
            # extend
            cov[-1][1] = sf[i]._end
        else:
            # start new
            cov.append([sf[i]._start, sf[i]._end])
    return cov

##INTERNALS##

known_attributes =set(['ID','Name','Alias','Parent','Target','Gap','Derives_from','Note','Dbxref','Ontology_term'])

def printWarning(string,fh=sys.stderr):
    fh.write('\n## WARNING ## ' + string + '\n')
    return

def suicide(string):
    sys.exit('\n## FATAL ERROR ## ' + string + '\n')
    return

def strandChar(st):
    if st in set(['+','-','.']):
        return st
    try:
        assert isinstance(st, int)
        if st > 0:
            return '+'
        elif st < 0:
            return '-'
        else:
            return '.'
    except:
        raise

def strandInt(st):
    if isinstance(st,int):
        return st
    try:
        assert st in set(['+','-','.'])
        if st == '+':
            return 1
        elif st == '-':
            return -1
        else:
            return 0
    except:
        raise

def binSearchRange(A,searchtuple,maxl):
    lowerIndex = None
    higherIndex = None

    loTup = (searchtuple[0], searchtuple[1] - maxl) # chr,start
    # lower boundary
    imin = int(0)
    imax = len(A)
    while imax > imin:
        imid = (imin + imax) / int(2)
        if (A[imid]._seqid, A[imid]._end) < loTup:
            imin = imid + 1
        else:
            imax = imid
    if imax == imin:
        lowerIndex = imin
    else:
        lowerIndex = 0

    # upper boundary
    hiTup = (searchtuple[0],searchtuple[2]) # chr,end
    imin = int(0)
    imax = len(A)
    while imax > imin:
        imid = (imin + imax) / int(2)
        if (A[imid]._seqid, A[imid]._start) < hiTup:
            imin = imid + 1
        else:
            imax = imid
    if imax == imin:
        higherIndex = imin
    else:
        higherIndex = len(A)

    # return range
    return (lowerIndex, higherIndex)

##CLASSES##

class Generic:
    instanceCount = 0
    @classmethod
    def incrementCount(self):
        self.instanceCount += 1
    
    def __init__(self,l=''):
        f = l.split('\t')
        try:
            f[8]
        except:
            print '## STRANGE LINE ##', l
            raise
        self._seqid = f[0]
        self._source = f[1]
        self._type = f[2]
        self._start = int(f[3]) if f[3] else None
        self._end = int(f[4]) if f[4] else None
        self._score = f[5]
        self._strand = f[6]
        self._phase = f[7]
        self._attributes = {}
        attrib = f[8].split(';')
        for a in attrib:
            if len(a) > 1:
                self.setAttribute(a.rstrip())
        # sortorder by type
        sortorder = { 'gene':1, 'mRNA':2, 'exon':3, 'three_prime_utr':4, 'start_codon':5, 'CDS':6, 'stop_codon':8, 'three_prime_utr':9 } # 7 os other
        self._sortorder = sortorder[self._type] if sortorder.has_key(self._type) else 7
        # invisible stuff
        self._flags = {}
        self._parent = None
        self._childs = [ ]
        # unique numbering
        self._number = self.instanceCount
        self.incrementCount()

    def __lt__(self,e):
        '''for sorting'''
        # LT
        if self._seqid < e._seqid:
            return True
        elif self._seqid == e._seqid:
            if self._start < e._start:
                return True
            elif self._start == e._start:
                if self._sortorder < e._sortorder:
                    return True
                elif self._sortorder == e._sortorder:
                    if self._end < e._end:
                        return True
                    else:
                        return False
        return False

    # a tuple that can be used for binary sorting
    def forSorting(self):
        return (self._seqid,self._start,self._end)

    def __str__(self):
        return "<Generic annotation chromosome:%s start:%d end:%d strand:%s>" % (self._seqid, self._start, self._end, self._strand)

    def __repr__(self,f=''):
        att = []
        for a in sorted(self._attributes.keys()):
            att.append('='.join([a,self._attributes[a]]))
        flg = []
        for f in self._flags.keys():
            if self._flags[f]: flg.append(f)
        if len(flg) > 0:
            att.append('FLAGS=' + ','.join(flg))
        return '\t'.join([ self._seqid, self._source, self._type, str(self._start), str(self._end), str(self._score), self._strand, self._phase]) + '\t' + ';'.join(att)

    def gtfstring(self):
        att = []
        for a in sorted(self._attributes.keys()):
            att.append(' '.join([ a, '"' + str(self._attributes[a]) + '"']))
        return '\t'.join([ self._seqid, self._source, self._type, str(self._start), str(self._end), self._score, self._strand, self._phase]) + '\t' + '; '.join(att)

    def ensgtfstring(self):
        ordering = { 'gene_id':1,'transcript_id':2,'exon_number':3,'protein_id':4 }
        att = []
        for a in sorted(self._attributes.keys()):
            if a in set(ordering.keys()):
                att.append(' '.join([ a, '"' + str(self._attributes[a]) + '"']))
        att = sorted(att,key=lambda x: ordering[x[:x.find(' ')]])
        return '\t'.join([ self._seqid, self._source, self._type, str(self._start), str(self._end), self._score, self._strand, self._phase]) + '\t' + '; '.join(att) + ';'

    def setUniqueIdent(self,prefix=''):
        if self._attributes.has_key('ID') and self._attributes['ID']:
            printWarning('Instance ID will be overwritten!')
        self.setAttribute('ID=' + str(prefix) + str(self._number))
        return self.getAttribute('ID')
    ## GETSET
    #attributes
    def setAttribute(self,a,warn=False):
        e = a.split('=')
        self._attributes[e[0]] = e[1]
        if warn and e[0] not in  known_attributes:
            printWarning('unknown attribute (' + e[0] + ')')
        return

    def hasAttribute(self,a,isSet=True):
        if isSet:
            return self._attributes.has_key(a) and self._attributes[a]
        else:
            return self._attributes.has_key(a)
    def getAttribute(self,a):
        try:
            return self._attributes[a]
        except:
            printWarning('attribute (' + a + ') does not exist')
            raise
            return

    def getAttributeString(self):
        att = []
        for a in self._attributes.keys():
            att.append('='.join([a,self._attributes[a]]))
        return ';'.join(att)

    def removeAttributesBut(self,keep=[]):
        newAttributes = {}
        for a in range(0,len(keep)):
            if self._attributes.has_key(keep[a]):
                newAttributes[keep[a]] = self._attributes[keep[a]]
        self._attributes = newAttributes
        return
    def length(self):
        return self._end - self.start(0)

    def start(self,base=1):
        if base == 0:
            return self._start - 1
        else:
            return self._start

    def end(self):
        return self._end

    def offset(self,offset):
        self._start += offset
        self._end += offset
        return self
    #flags
    def setFlag(self,flag):
        self._flags[flag] = True
        return self._flags[flag]

    def unsetFlag(self,flag):
        if self._flags.has_key(flag):
            self._flags[flag] = False
        else:
            printWarning('flag (' + flag + ') cannot be unset')
            return
        return self._flags[flag]

    def hasFlag(self,flag):
        if self._flags.has_key(flag):
            return self._flags[flag]
        else:
            return False
    #strand
    def getStrand(self):
        return self._strand
    def getIntStrand(self):
        return strandInt(self._strand)
    def setStrand(self,st):
        self._strand = strandChar(st)
        return
    #child/parent
    def createParent(self,parent_ids,annotationType,debug=False):
        # carboncopy
        line = repr(self)
        parent = Generic(line)
        # remove everything but parent_id
        parent.removeAttributesBut(parent_ids)
        # set id
        parent.setAttribute('ID=' + parent.getAttribute(parent_ids[0]) + '.' + self._seqid)
        # set type
        parent._type = annotationType
        # remove phase
        parent._phase = '.'
        # return newly created parent
        return parent
        
    def setParent(self,p):
        assert isinstance(p, Generic)
        if self._attributes.has_key('Parent') and self._parent:
            if self._parent is not p:
                printWarning('parent changed')
        self._parent = p
        self._attributes['Parent'] = self._parent.getAttribute('ID')
        return

    def hasParent(self):
        return self._attributes.has_key('ID')

    def addChild(self,e,check=True):
        assert isinstance(e, Generic)
        self._childs.append(e)
        self._childs[-1].setParent(self)
        if check:
            # update values
            if not self._seqid and e._seqid:
                self._seqid = e._seqid
            else:
                try:
                    assert self._seqid == e._seqid
                except:
                    self.setFlag('chimeric')
             
            if not self.hasFlag('chimeric'):
                if not self._start or e._start < self._start:
                    self._start = e._start
                if not self._end or e._end > self._end:
                    self._end = e._end
        return

    def getGrandChilds(self):
        gc = []
        for c in self._childs:
            gc.extend(c._childs)
        # make it unique and return
        return list(set(gc))
    
    ## BASIC STUFF
    def overlap(self,e):
        if e._seqid == self._seqid:
            if e._end < self._start or self._end < e._start:
                return False
            else:
                return True
        else:
            return False

    def joint(self,e):
        if e._seqid == self._seqid:
            if e._end < self._start - 1 or self._end + 1 < e._start:
                return False
            else:
                return True
        else:
            return False

    def extend(self,e):
        try:
            assert self.overlap(e)
        except:
            sys.stderr.write('Unable to extend non-overlapping exons')
            raise
        # adapt field values
        if self._type != e._type:
            self._type = ':'.join(list(set(self._type.split(':') + e._type.split(':'))))
        if self._source != e._source:
            self._source = 'various'
        self._score = '.'
        if self._strand != e._strand:
            self._strand = '.'
        # adapt phase
        if e._start < self._start and self._phase != '.':
            self._phase = (self._phase + self._start - e._start) % 3
        # reset uncommon attributes
        for a in self._attributes.keys():
            if not e._attributes.has_key(a) or e._attributes[a] != self._attributes[a]:
                del self._attributes[a]
        # remove all flags
        self._flags = {}
        # extend coordinates
        if e._end >= self._end:
            self._end = e._end
        if e._start <= self._start:
            self._start = e._start
        return

    ## ensembl filtering options (for eliminating out of phase projections)
    def cdsExons(self):
        """returns coding exons"""
        cds = []
        for c in self._childs:
            if c._type == 'CDS':
                cds.append(c)
        return cds

    def phaseBreak(self):
        """returns last cdsexons offset with correct phase"""
        cds = self.cdsExons()
        if self._strand == '-':
            cds.sort()
            cds.reverse()
        elif self._strand == '+':
            cds.sort()
        else:
            printWarning("phasecheck attempted without knowning the transcripts orientation")
        phase = 0
        cdsnumber = -1
        for c in cds:
            cdsnumber += 1
            if phase != int(c._phase):
                return cdsnumber
            phase = (int(c._phase) + ((3 - c.length()) % 3)) % 3
        if phase != 0:
            return cdsnumber
        return cdsnumber + 1

    def checkPhase(self):
        """checks if CDS phase and order is consistent"""
        cds = []
        # get CDS childs
        for c in self._childs:
            if c._type == 'CDS':
                cds.append(c)
        # sort them
        if self._strand == '-':
            cds.sort()
            cds.reverse()
        elif self._strand == '+':
            cds.sort()
        else:
            printWarning("phasecheck attempted without knowning the transcripts orientation")
        # check Phase
        ##resu = ''
        phase = 0
        for c in cds:
            if phase != int(c._phase):
                ##resu += '0'
                return False
            else:
                ##resu += '1'
                pass
            phase = (int(c._phase) + ((3 - c.length()) % 3)) % 3
        # phase at the end (should be zero)
        if phase != 0:
            ##resu += '0'
            return False
        else:
            ##resu += '1'
            pass
        ##if '0' in resu: return False
        ##else: return True
        return True

    def addSubId(self,sub,ids=['ID'],sep='$'):
        for i in ids:
            if self._attributes.has_key(i):
                self._attributes[i] += sep + sub
            else:
                printWarning('Cannot add subid to %s in %s' % (i, str(self)))
        return self
    def strictOverlap(self,e,req):
        # consider strandedness when overlapping (and sureStrand)
        if self.overlap(e) and abs(self.sureStrand(req) - e.sureStrand(req)) < 2:
            return True
        return False

    def strictJoint(self,e,req):
        # consider strandedness when overlapping (and sureStrand)
        if self.joint(e) and abs(self.sureStrand(req) - e.sureStrand(req)) < 2:
            return True
        return False

    def sureStrand(self,req='spliced'):
        if self.hasFlag(req):
            if self._strand == '+':
                return 1
            elif self._strand == '-':
                return -1
            else:
                return 0
        else:
            return 0

    def majStrand(self):
        s = {}
        for c in self._childs:
            if not s.has_key(c._strand):
                s[c.getStrand()] = 0
            s[c.getStrand()] += 1
        ss = sorted(s, key=s.get, reverse=True)
        return ss[0]

    def implicitStrand(self,required=None):
        st = 0
        for c in self._childs:
            try:
                assert abs(st - c.sureStrand(required)) < 2
            except:
                printWarning('Cannot get implicit strand, childs are ambiguous')
                self._strand = '.'
                return
            if c.sureStrand(required) != 0 and st == 0:
                st = c.sureStrand(required)
        return st

    def flattenChilds(self,joint=False):
        flatChilds = []
        for e in sorted(self._childs):
            if len(flatChilds) == 0:
                flatChilds.append(e)
            else:
                if flatChilds[-1].overlap(e) or (joint and flatChilds[-1].joint(e)):
                    flatChilds[-1].extend(e)
                    if self._strand != '.':
                        flatChilds[-1].setStrand(self._strand)
                else:
                    flatChilds.append(e)
        # reassign IDs
        exonId = 1
        for e in sorted(flatChilds):
            e.setAttribute('ID=' + e.getAttribute('Parent') + '.' + str(exonId))
            exonId += 1
        # set new childs
        self._childs = flatChilds
        return len(self._childs)

    def write(self,handle,childs=False,parent=True,depth=1,formatting='GFF'):
        if parent:
            if formatting=='GTF':
                print >> handle, self.gtfstring()
            else:
                print >> handle, repr(self)
        if childs:
            #print self.getAttribute('ID'), 'writing childs', len(self._childs)
            for c in sorted(self._childs):
                if depth > 1:
                    c.write(handle, True, True, depth - 1,formatting)
                elif formatting=='GTF':
                    print >> handle, self.gtfstring()
                else:
                    print >> handle, repr(self)
        return

    def bottomLevelOverlap(self,other):
        if self._end < other._start or self._start > other._end or self._strand != other._strand:
            return 0
        # overlap lists (ordered so can break when needed)
        c = sorted(covered(self.bottomLevel()) + covered(other.bottomLevel()))
        covScore = 0
        for i in range(1,len(c)):
            if c[i]._start <= c[i-1]._end:
                covScore += c[i-1]._end - c[i]._start + 1
        return covScore

    def bottomLevel(self):
        # returns an list of tuples of stretches covered by bottomlevel
        if len(self._childs) > 0:
            thisChilds = self._childs
            while len(thisChilds) > 0:
                nextChilds = []
                for c in thisChilds:
                    nextChilds.extend(c._childs)
                if len(nextChilds) == 0:
                    return thisChilds
                else:
                    thisChilds = nextChilds
            assert False
        return self

class GTF(Generic):
    def __init__(self,l,fields=set(['gene_id'])):
        f = l.split('\t')
        g = f[:8]
        att = f[8].split(';')
        attributes = []
        ID = []
        for a in att:
            if len(a) > 3:
                a = a.strip() # strip whitspace
                aSplit = a.split(' "') # split key and value
                aSplit[1] = aSplit[1][:-1] # strip quotes
                attributes.append('='.join(aSplit))
                if aSplit[0] in fields:
                    ID.append(aSplit[1])
        # put an ID
        attributes.append('ID=' + '|'.join(ID))
        g.append(';'.join(attributes)) # g[8]
        GFFline = '\t'.join(g)
        Generic.__init__(self,GFFline)

class Gene(Generic):
    def __init__(self,f):
        f[0] = f[0] if f[0] else ''
        f[1] = f[1] if f[1] else 'unknown'
        f[2] = f[2] if f[2] else 'unknown'
        f[3] = f[3] if f[3] else ''
        f[4] = f[4] if f[4] else ''
        f[5] = f[5] if f[5] else '.'
        f[6] = f[6] if f[6] else '.'
        f[7] = f[7] if f[7] else '.'
        attrib = []
        for k in f[8].keys():
            attrib.append(k + '=' + f[8][k])
        f[8] = ';'.join(attrib)
        Generic.__init__(self,'\t'.join(f))

class Parent(Generic):
    def __init__(self,array,ptype):
        # get boundaries
        boundaries = [None,None]
        chromosome = None
        strand = None
        strandcount = [0,0]
        attribs = copy.copy(array[0]._attributes) # dic (shallow copy sufficient)
        source = array[0]._source
        ambiguousAttributes = set(['Parent']) # parent wont obtain its own parent reference
        for a in array:
            if not boundaries[0] or a.start(1) < boundaries[0]:
                boundaries[0] = a.start(1)
            if not boundaries[1] or a.end > boundaries[1]:
                boundaries[1] = a.end()
            if chromosome and chromosome != a._seqid:
                suicide('Cannot build parent for multiple seqIDs')
            chromosome = a._seqid
            if strand and strand != a._strand:
                suicide('Cannot build parent for ambiguous strands')
            strand = a._strand
            if source and source != a._source:
                printWarning('Parent will be built from mixed sources')
            source = a._source
            # determine ambiguous attributes
            for att in a._attributes.keys():
                if not attribs.has_key(att) or (attribs.has_key(att) and attribs[att] != a._attributes[att]):
                    ambiguousAttributes.add(att)
        # remove ambiguous attributes
        for a in ambiguousAttributes:
            if attribs.has_key(a):
                del attribs[a]
        # build entry
        new = [None]*9
        attrib = []
        new[0] = chromosome
        new[1] = source
        new[2] = ptype
        new[3] = str(boundaries[0])
        new[4] = str(boundaries[1])
        new[5] = '.'
        new[6] = strand
        new[7] = '.'
        for k in attribs.keys():
            attrib.append(k + '=' + attribs[k])
        new[8] = ';'.join(attrib)
        Generic.__init__(self,'\t'.join(new))

class Projection(Generic):
    def __init__(self,origin,targetCoordinates,targetSpecies):
        t = None
        if targetCoordinates[1] <= targetCoordinates[2]:
            # positive
            t = (targetCoordinates[0],targetCoordinates[1],targetCoordinates[2],targetCoordinates[3])
        else:
            # negative
            t = (targetCoordinates[0],targetCoordinates[2],targetCoordinates[1],targetCoordinates[3])
        f = [ str(t[0]), origin._source, targetSpecies + '_PROJECTION', str(t[1]), str(t[2]), '.', str(strandChar(t[3])), '.', origin.getAttributeString() ]
        Generic.__init__(self,'\t'.join(f))

class Fragment(Generic):
    def __init__(self,a,b,c,d):
        gf = []
        gf.append(str(a))
        gf.append('unknown')
        gf.append('fragment')
        gf.append(str(b))
        if c == None or c == 'None': gf.append(str(9999999999))
        else: gf.append(str(c))
        gf.append('.')
        gf.append('.')
        gf.append('.')
        gf.append('ID=' + str(d) + ';')
        Generic.__init__(self,'\t'.join(gf))
