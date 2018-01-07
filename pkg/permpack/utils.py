import sys
        
def normalize(l, n):

        """ make ordered numbers into a permutation (not instance of  Permutation)
        
        INPUT
        l: list of distinct integers
        n: length of lst
        
        EXAMPLE
        
        sage: normalize([2,5,3],3)
        sage: [1,3,2]
        """
        lst = range(n)
        lst.sort(key=lambda i: l[i])
        normal_lst = range(n)
        for i, x in enumerate(lst):
            normal_lst[x] = i + 1
        return normal_lst
"""
def complement(pos,m):

        return list(set(range(m))-set(pos))


def normal_subperm(perm, tp, m, k, n):
        r"
        perm: to be permutation, now only m distinct integers
        tp: distinguished positions
        n: length of original perm
        k: length of tp
        m: length of perm
        "
        
        MAXVAL = 999999999999 # hopefully big enough to be MAXVAL
        normal_perm = [0 for x in perm]
        normal_tp = [0 for x in tp]
        tset = Set(tp)

        for val in range(1,m+1):
                imin = perm.index(min(perm))
                normal_perm[imin] = val
                if imin in tset:
"""                     
