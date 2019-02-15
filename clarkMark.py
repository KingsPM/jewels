#!/usr/bin/env python

import time

def clarkMark():
    start_time = time.time()
    sum = 0 
    for i in range(1000000):
	sum = sum + i
    timeElapsed =  time.time() - start_time
    return 1.0 / timeElapsed

if __name__=="__main__":
    cm = []
    for i in range(10):
        cm.append(clarkMark())
        print i, "ClarkMark", cm[-1]
        print i, "PiPower", cm[-1]/0.586
        
