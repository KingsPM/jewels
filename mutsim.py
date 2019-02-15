#!/usr/bin/env python

'''
mutsim.py
Ver primitive simulation for neutral sequence evolution over time
'''

import sys
import string
import random
from copy import deepcopy


def randomseq(size=1000, gccontent=0.5, at=['A', 'T'], gc=['G', 'C']):
    gclen = int(size*gccontent)
    atlen = size - gclen
    bases = [random.choice(gc) for x in range(gclen)] + [random.choice(at) for x in range(atlen)]
    random.shuffle(bases)
    return bases


def identity(a, b):
    ident = 0
    for i in range(len(a)):
        if a[i] == b[i]:
            ident += 1
    return ident/float(len(a))


def samplepositions(l, n):
    return [random.randint(0, l-1) for i in range(n)]


def mutate(x, k=1.0/3.0):
    # transitions are 2x more likely
    if k < random.random():
        r = random.random()
        # tranversion
        if x == 'A':
            return 'C' if r < 0.5 else 'T'
        elif x == 'G':
            return 'C' if r < 0.5 else 'T'
        elif x == 'C':
            return 'A' if r < 0.5 else 'G'
        elif x == 'T':
            return 'A' if r < 0.5 else 'G'
    else:
        # transition
        if x == 'A':
            return 'G'
        elif x == 'G':
            return 'A'
        elif x == 'C':
            return 'T'
        elif x == 'T':
            return 'C'
    return '?'

if __name__ == "__main__":

    seqlen = 1000000

    rate = 1e-9 if len(sys.argv) < 2 else float(sys.argv[1])
    step = 1000000
    evotime = 200

    final = 0.3 if len(sys.argv) < 3 else float(sys.argv[2])

    initial = randomseq(seqlen)
    mutable = deepcopy(initial)
    mutations = int(rate*float(step)*float(seqlen))
    print mutations
    for i in range(evotime):
        currenttime = i*step
        # get mutable positions
        mutpos = samplepositions(seqlen, mutations)
        for p in mutpos:
             mutable[p] = mutate(mutable[p])
        # get identity
        ident = identity(initial,mutable)
        print i*step, ident
        if 1.0-ident > final:
            break


