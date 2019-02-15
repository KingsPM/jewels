#!/usr/bin/env python

def binsearch(key, lst):
    try: mid = lst[len(lst)/2]
    except: return []
    else:
        ov = []
        if overlaps(key,mid):
            ov.append(mid)
            for k in lst[len(lst)/2+1:]:
                if overlaps(key,k):
                    ov.append(k)
                else:
                    break
            for i in range(len(lst)/2-1,-1, -1):
                if overlaps(key,lst[i]):
                    ov.append(lst[i])
                else:
                    break
            return ov
        else:
            return binsearch(key,lst[:len(lst)/2])+binsearch(key, lst[len(lst)/2+1:])

def overlaps(A,B):
    if A[1] < B[0]: return False
    elif A[0] > B[1]: return False
    else: return True


if __name__ == "__main__":

    test = [(i,i+10) for i in range(10,50,3)]

    key = (20,40)
    
    print 'query', key
    print 'test ', test
    print sorted(binsearch(key,test))
