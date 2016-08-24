from sage.combinat.permutation import Permutation, Permutations, PatternAvoider
from sage.combinat.combination import Combinations
from sage.arith.all import factorial
from sage.rings.all import Integer
from utils import *

class PermFlag:
    
    def __init__(self, perm, tpos=list()):
        """
        Returns a PermFlag object.

        INPUT:
        
        - perm: permutation as a list of int and letters or string of those
        - tpos: type positions (list of ints).
        
        EXAMPLE:
        
        sage: pf1 = PermFlag([2,3,4,1,5], [1,2,4])
        sage: pf2 = PermFlag("123456789bca")
        """
        # if perm is a string
        if isinstance(perm, basestring):
            permlist = list(perm)
            intlib = list("123456789")
            charlib = list("abcdefghijklmnopqrstuvxyz")
            perm_candidate = []
            for x in perm:
                if x in intlib:
                    perm_candidate.append(int(x))
                elif x in charlib:
                    i = charlib.index(x)
                    perm_candidate.append(10+i)
                else:
                    raise ValueError("Enter a proper permutation!")
            self.perm = Permutation(perm_candidate)
            
        # if perm is integer list
        elif all(isinstance(x, Integer) or isinstance(x,int) for x in perm):
            self.perm = Permutation(perm)
        else:
            raise ValueError("Pattern needs to be either string of integers or list of integers.")
        
        self.N = len(perm)
        self.tp = self._induced_subperm(tpos)
        self.tp_pos = tpos
        self.tn = len(tpos)

        
    def __eq__(self, pf):
        """ Override '==' operator on PermFlag class."""
        if isinstance(pf,self.__class__):
            return (pf.perm == self.perm) and (self.tp_pos == pf.tp_pos)
        else:
            return False


    def __repr__(self):

        if self.tn > 0:
            return self.perm.__str__()+self.tp_pos.__str__()

        return self.perm.__str__()    

    def _to_str(self):

        return "".join([str(x) for x in self.perm])

    def _induced_subperm(self, positions=list()):

        """ Does NOT check whether the whole type is in subperm. """

        return Permutation(normalize([self.perm[i] for i in positions], len(positions)))


    def _induced_subpattern(self, positions=list()):

        """Return subpattern on given positions. Normalized."""

        return normalize([self.perm[i] for i in positions], len(positions))

    def subperm_density(self, perm):
        """
        Return density of perm in self.

        INPUT:
        
        - perm:    permutation of type PermFlag

        EXAMPLE:

        sage: P = PermFlag("1234765")
        sage: P.subperm_density(PermFlag("123"))
        22/35
        """

        combs = Combinations(self.N, perm.N)
        counter = 0
        for c in combs:
            if self._induced_subpattern(c) == perm.perm:
                counter += 1
        return Integer(counter)/combs.cardinality()
