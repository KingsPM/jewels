#1/usr/bin/env python

'''
parser for newick tree files
builds directed graph of tree (leaf->root)
'''

import sys
import time
import re
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from collections import Counter

# parenthesis matcher
parma = re.compile(r"\(([^\()]+)\)")

class Newick:
    def __init__(self,treeline):
        self.root = None
        self.network = self._parse(treeline.rstrip())
        return

    def _parse(self,nw):
        # sanity check
        a = nw.count('(')
        try:
            assert a == nw.count(')') and nw[-1:] == ';'
        except:
            print >> sys.stderr, "## ERROR ## invalid newick tree"
            print >> sys.stderr, "## TREE ###", nw, "<"
            raise
        # parse newick tree structure
        network = nx.DiGraph()
        ancNumber = 1
        atRoot = False
        m = parma.search(nw)
        while m and not atRoot:
            # rebuild string
            s =  m.span()
            anc = 'node'+str(ancNumber)
            if nw[:s[0]].count('(') == 0:
                atRoot = True
                self.root = anc
            else:
                # not root
                nw = nw[:s[0]] + anc + nw[s[1]:]
            # update graph
            descendants = m.group(1).split(',')
            network.add_node(anc)
            for d in descendants:
                network.add_edge(anc,d)
            # rematch string
            m = parma.search(nw)
            ancNumber += 1
        return network

    def __repr__(self):
        return self.newickString('')  # return newick string without labels

    def newickString(self,labels='label',verbose=False):
        '''resolves graph structure into tree'''
        G = self.network
        visited = set()
        # start from root
        dd = sorted(G[self.root].keys())
        visited.add(self.root)
        nw = '(' + ','.join(dd) + ')'
        descendants = [ n for n in dd if len(G.neighbors(n)) > 1]
        # get labels
        nodelabels = nx.get_node_attributes(self.network,labels)
        # solve all non-terminal descendants
        while descendants:
            # iterate descendants
            d = descendants.pop(-1)
            dd = list(set(G[d].keys()) - visited)
            descendants += [ n for n in dd if len(G.neighbors(n)) > 1 ]
            # manipulate newick string
            substitution = '(' + ','.join(dd) + ')'
            # add label
            try:
                nodelabels[d]
            except:
                pass
            else:
                substitution += ' ' + nodelabels[d]
                del nodelabels[d]
            # update newick string
            if verbose:
                print >> sys.stderr, ">>", nw, d, '=>', substitution
            matchstring = r'\b'+d+r'\b'
            try:
                assert len(re.findall(matchstring,nw)) == 1
            except:
                print >> sys.stderr, '## ERROR ## cannot substitute node %s with %s' % (d, substitution)
                raise
            else:
                nw = re.sub(matchstring, substitution, nw, count=1)
                if verbose:
                    print >> sys.stderr, "--", nw
                visited.add(d)

        # add labels to terminal branches
        if nodelabels:
            for tip,lab in nodelabels.items():
                if verbose:
                    print >> sys.stderr, ">>", nw, tip, 'labeled', lab
                matchstring = r'\b'+tip+r'\b'
                substitution = tip+' '+lab
                try:
                    assert len(re.findall(matchstring,nw)) == 1
                except:
                    print >> sys.stderr, '## ERROR ## cannot add label %s to tip %s' % (lab, tip)
                    raise
                else:
                    nw = re.sub(matchstring, substitution, nw, count=1)
                    if verbose:
                        print >> sys.stderr, "--", nw
        return nw + ';'

    def _checkNode(self,n):
        try:
            self.network[n]
        except:
            print >> sys.stderr, "## ERROR ## Node %s not found" % (n)
            raise
        return


    def viz(self,filename='tree.png'):
        #pos=nx.spring_layout(self.network,scale=1) #default to scale=1
        pos=nx.spectral_layout(self.network,scale=1) #default to scale=1
        nx.draw(self.network,pos)
        plt.savefig(filename)
        plt.clf()
        return

    def allLeaves(self):
        return [ n for n in self.network.nodes() if self.network.degree(n) == 1 ]

    def mrca(self,*leaves):
        # sanity check
        try:
            assert set(leaves).issubset(set(self.allLeaves()))
        except:
            print >> sys.stderr, "## ERROR ## Check leaves for mrca"
            raise
        # get shortest paths to root
        hitnodes = []
        for l in leaves:
            hitnodes += nx.shortest_path(self.network, self.root, l)[:-1]
        counted = Counter(hitnodes)        
        countedmax = [node for node, count in counted.items() if count == len(leaves)]
        return str(min(countedmax))

    def leaves(self,node):
        try:
            self.network[node]
        except:
            print >> sys.stderr, "## ERROR ## there is no node named", node, "in the tree"
            raise
        return [ l for l in self.allLeaves() if nx.has_path(self.network, node, l) ]

    def isMonophyletic(self,*leaves):
        commonAncestor = self.mrca(*leaves)
        commonAncestorLeaves = set(self.leaves(commonAncestor))
        if len(set(leaves).symmetric_difference(commonAncestorLeaves)) == 0:
            return True
        return False

    def unroot(self,verbose=False):
        '''unroots tree to comply with PAML standards'''
        oldroot = self.root
        childs = self.network[self.root].keys()
        if len(childs) > 2:
            # already unrooted
            pass
        elif len(childs) == 2:
            for c in childs:
                if self.network.degree(c) == 1:
                    pass
                else:
                    self.network.remove_node(self.root)
                    self.root = c
                    i = 0 if self.root is childs[1] else 1  ## decides direction of edge (root->leaves)
                    self.network.add_edge(self.root,childs[i])
                    break
            if verbose:
                print >> sys.stderr, "UNROOTING SUCCESS: changed headnode", oldroot, '->', self.root
        return self

    def branchLabel(self,n,value):
        nx.set_node_attributes(self.network, 'label', {n: value})
        ancestor = self.ancestor(n)
        nx.set_edge_attributes(self.network, 'label', {(ancestor,n): value})
        return

    def ancestor(self,n):
        self._checkNode(n)
        # get ancestor
        try:
            assert len(self.network.reverse()[n].keys()) == 1
        except:
            print >> sys.stderr, "## ERROR ## Node %s has %d ancestors (should be exactly one)" % (n,len(self.network.reverse()[n].keys()))
            raise
        # return ancestor node
        return self.network.reverse(copy=True)[n].keys()[0]

    def descendants(self,n):
        self._checkNode(n)
        return self.network[n].keys()

    # reroot (changes edge directions, removes root node and reconnects childs, inserts newroot node)


if __name__=="__main__":
    sampletrees = [ '(((A,B),C),((D,E),(F,G,(I,(J,K)))));', '(((A,(B,K)),C),((((D,E),(F,G)),I),J));' ]
    for sampletree in sampletrees:
        tree = Newick(sampletree)
        print 'INPUT', sampletree
        print '\tPARSED', repr(tree)
        print '\tLEAVES', tree.allLeaves()
        print '\tMRCA(ABC)', tree.mrca('A','B','C')
        print '\tTIPS(node5)', tree.leaves('node5')
        print '\tMONOPHYLETIC(AKB)', tree.isMonophyletic('A','K','B')
        print '\tUNROOTED', repr(tree.unroot())
        #tree.viz('sample.png')

    if tree.isMonophyletic('A','K','B'):
        mrca = tree.mrca('A','K','B')
        labelnumber = 1
        tree.branchLabel(mrca, '$'+str(labelnumber))
        descendants = tree.descendants(mrca)   
        for d in descendants:
            labelnumber += 1
            tree.branchLabel(d, '#'+str(labelnumber))
        print "withLabels", tree.newickString()

    ### tree ordering must be like sequence ordering
