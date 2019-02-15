#!/usr/bin/env python

'''
module with simple machine learning functions:
DBSCAN clustering
'''

from math import sqrt, pow
import sys


class DBSCAN:
    # Density-Based Spatial Clustering of Application with Noise -> http://en.wikipedia.org/wiki/DBSCAN
    def __init__(self, esp=4, MinPts=2):
        self.name = 'DBSCAN'
        self.DB = []  # Database
        self.esp = esp  # neighborhood distance for search
        self.MinPts = MinPts  # minimum number of points required to form a cluster
        self.cluster = []

    def DBSCAN(self):
        for P in self.DB:
            if not P.visited:
                #for each unvisited point P in dataset
                P.visited = True
                NeighborPts = self.regionQuery(P)
                if len(NeighborPts) < (self.MinPts - 1):
                    #that point is a noise (no neighbors)
                    P.isnoise = True
                    #print P.show(), 'is a noise'
                else:
                    self.cluster.append([])
                    self.expandCluster(P, NeighborPts)

    def expandCluster(self, P, neighbor_points):
        self.cluster[-1].append(P)
        P.membership = True
        for np in neighbor_points:
            if not np.visited:
                #for each point P' in NeighborPts
                np.visited = True
                NeighborPts_ = self.regionQuery(np)
                if len(NeighborPts_) >= self.MinPts:
                    for np_ in NeighborPts_:
                        neighbor_points.append(np_)
            if not np.membership:
                #if P' is not yet member of any cluster
                self.cluster[-1].append(np)
                np.membership = True
        return

    def regionQuery(self, P):
        #return all points within P's eps-neighborhood, except itself
        pointInRegion = []
        for p_tmp in self.DB:
            if P.dist(p_tmp) < self.esp and P is not p_tmp:
                pointInRegion.append(p_tmp)
        return pointInRegion


class Point:
    def __init__(self, coord, name=None, visited=False, isnoise=False, membership=False):
        self.coord = coord
        self.dim = len(coord)
        self.name = name
        self.visited = visited
        self.isnoise = isnoise
        self.membership = membership

    def dist(self, other):
        try:
            assert self.dim == other.dim
        except:
            print >> sys.stderr, "ERROR: incompatible dimensions %d != %d" % (self.dim, other.dim)
            raise
        return sqrt(sum([pow(other.coord[i]-self.coord[i], 2) for i in range(self.dim)]))

    def show(self):
        if self.name:
            return self.name
        return self.coord


if __name__=='__main__':
    #this is a mocking data just for test
    vecPoint = [
        Point((11, 3, 2)),
        Point((10, 4, 2)),
        Point((11, 5, 4)),
        Point((12, 4, 6)),
        Point((13, 5, 20)),
        Point((12, 6, 20)),
        Point((13, 5, 19)),
        Point((6, 10, 12)),
        Point((8, 10, 129)),
        Point((5, 12, 13)),
        Point((7, 12, 5)),
        Point((5, 5, 5))
    ]

    #Create object
    dbScan = DBSCAN()
    #Load data into object
    dbScan.DB = vecPoint
    #Do clustering
    dbScan.DBSCAN()
    #Show result cluster
    for i in range(len(dbScan.cluster)):
        print 'Cluster: ', i
        for j in range(len(dbScan.cluster[i])):
            print dbScan.cluster[i][j].show()
