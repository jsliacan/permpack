import os, sys, time, multiprocessing

from numpy import array

from sage.combinat.permutation import Permutation, Permutations, PatternAvoider
from sage.combinat.combination import Combinations
from sage.combinat.subset import Subsets, Set
from sage.arith.all import factorial, binomial
from sage.rings.all import Integer, Rational
from sage.doctest.util import Timer
from sage.misc.functional import numerical_approx
from sage.matrix.all import matrix

from perm_flag import *
from utils import *

NUMCPU = multiprocessing.cpu_count() # global var, why not

class PermProblem:

    

    def __init__(self, N, forbid=None, density_pattern=None):
        """
        Return an instance of a PermProblem.

        INPUT:

        - N:                         order of perms in admissible family
        - forbid:                    list of permutations (strings or lists(of ints))
        - density_pattern:           either a single perm, or list of pairs (perm, coef)

        EXAMPLE:

        sage: p = PermProblem(3, density_pattern=[1,3,2])
        Generating admissible permutations: OK.
        Generating types: OK.
        Generating flags: OK.
        Expressing density pattern as a linear combination of admissible permutations: OK.
        ------------------------------------------------
        Generated:
        6 admissible permutations.
        1 types of order 1, with [4] flags.
        ------------------------------------------------
        Loading flag_products from file: OK.
        """

        self.N = N

        self.type_orders = list()
        for i in range(self.N):
            if (self.N-i)%2==0:
                self.type_orders.append(i)

        self.num_types = [0 for x in self.type_orders]

        
        if density_pattern == None:
            self.density_pattern = PermFlag([2,1])
        else:
        
            self.density_pattern = self._set_density(density_pattern)

        if forbid == None:
            self.forbidden_perms = list()
        else:
            self.forbidden_perms = [PermFlag(x) for x in forbid]

        # ADMISSIBLE PERMS
        sys.stdout.write("Generating admissible permutations: ")
        sys.stdout.flush()
        self.admissible_perms = self._generate_admissible_perms()
        sys.stdout.write("OK.\n")
        sys.stdout.flush()
        
        # TYPES
        sys.stdout.write("Generating types: ")
        sys.stdout.flush()
        self.types = self._generate_types()
        sys.stdout.write("OK.\n")

        self._info = len(self.types)

        # FLAGS
        sys.stdout.write("Generating flags: ")
        sys.stdout.flush()
        # generate flags
        self.flags = self._generate_flags()
        sys.stdout.write("OK.\n")

        # COMPUTING DENSITIES
        sys.stdout.write("Expressing density pattern as a linear combination of admissible permutations: ")
        sys.stdout.flush()
        self.densities = self._compute_densities()
        sys.stdout.write("OK.\n")
        sys.stdout.flush()

        self.assumptions = False
        self.assumption_densities = list()

        # PRINTING INFO
        sys.stdout.write("------------------------------------------------\n")
        sys.stdout.write("Generated:\n%d admissible permutations.\n" % len(self.admissible_perms))
        counter = 0
        until = 0
        for i in range(len(self.num_types)):
            until += self.num_types[i]
            sys.stdout.write("%d types of order %d, with " % (self.num_types[i], self.type_orders[i]))
            flag_counts = list()
            while counter < until:
                flag_counts.append(len(self.flags[counter]))
                counter += 1
            sys.stdout.write(str(flag_counts))
            sys.stdout.write(" flags.\n")
        sys.stdout.write("------------------------------------------------\n")
        sys.stdout.flush()

        # COMPUTE FLAG PRODUCTS
        self.flag_products = list()

        # timer = Timer().start()
        self._compute_flag_products()
        # print timer.stop()

        

    def __str__(self):

        return "A Permutation Problem class. Read the code to know more."



    def _set_density(self, dp):
        """Parse linear combination of densities."""

        density = None
        
        if isinstance(dp[0],tuple):
            density = [(PermFlag(x[0]), Rational(x[1])) for x in dp]
        else:
            density = [(PermFlag(dp),Rational(1))]

        return density
 
        
    def _is_admissible(self, perm):
        """Check if perm is admissible."""

        if not isinstance(perm, PermFlag):
            raise ValueError("Input must be of a valid PermFlag object.")
        
        order = perm.N
        verdict = True
        
        for fp in self.forbidden_perms:
            perm_contains_fp = False
            for c in Combinations(order, fp.N):
                induced_subperm = perm._induced_subpattern(c)
                if induced_subperm == fp.perm:
                    perm_contains_fp = True
                    break
            if perm_contains_fp:
                verdict = False
                break
            
        return verdict

    
    def _generate_admissible_perms(self, order=None):

        """Generate admissible permutations of order N (free from forbidden patterns)."""
        admissible_perms = list()

        if order == None:
            order = self.N

        for x in Permutations(order):
            candidate = PermFlag(x)
            if self._is_admissible(candidate):
                admissible_perms.append(candidate)

        return admissible_perms

    
    def _generate_types(self):
        
        types = list()
        type_orders = list()

        counter = 0
        for n in self.type_orders:
            types_n = self._generate_admissible_perms(n)
            self.num_types[counter] = len(types_n)
            types.extend(types_n)
            counter += 1
        
        return types
    
    
    def _generate_flags_slow(self):

        sys.stdout.write("Generating flags (SLOW VERSION): ")
        sys.stdout.flush()
        
        flags = list()

        sys.stdout.write("Generating flags: ")
        
        for t in self.types:
            t_flags = list()
            n = int((self.N-t.N)/2 +t.N) # must be casted as int
            allperms = [PermFlag(x) for x in Permutations(n)]
            S = Subsets(range(n),t.N)

            for p in allperms: # check each permutation 
                for s in S: # and all possible placings of type in it
                    subp = PermFlag(normalize([p.perm[i-1] for i in s], t.N))
                    if subp == t:
                        t_flags.append(PermFlag(p.perm, list(s)))

            flags.append(t_flags)

        sys.stdout.write("OK.\n")
        sys.stdout.flush()
        
        return flags
    

    def _generate_flags(self):

        """Return list of lists of admissible flags for each type (one list of flags for each type)."""
        flags = list()

        for t in self.types:
            t_flags = list() # list of flags on type t
            tperm = t.perm
            m = int((self.N+t.N)/2) # must be casted as int
            rm = range(m)
            srm = Set(rm)
            petals = Permutations(m-t.N)
            subsets = Subsets(range(m), t.N) # will serve as values and positions
            cosubsets = [srm-s for s in subsets]
            for petal in petals: # now type and petal are both fixed
                for s in range(len(subsets)):
                    positions = subsets[s]
                    copositions = cosubsets[s]
                    for v in range(len(subsets)):
                        values = subsets[v]
                        covalues = cosubsets[v]

                        newflag = range(m)
                        for i in range(t.N):
                            newflag[positions[i]] = values[tperm[i]-1]+1
                        for j in range(m-t.N):
                            newflag[copositions[j]] = covalues[petal[j]-1]+1

                        if self._is_admissible(PermFlag(newflag)):
                            t_flags.append(PermFlag(newflag,list(subsets[s])))
                            
            flags.append(t_flags)

        return flags



    def _compute_densities(self):
        
        densities = [0 for x in range(len(self.admissible_perms))]

        for densperm, coeff in self.density_pattern:

            if self.N < densperm.N: # too big to fit, then skip
                continue
            
            combs = Combinations(range(self.N), densperm.N)
            numcombs = combs.cardinality()
            
            # for each admissible permutation, compute the density of densperm in it
            for i in range(len(self.admissible_perms)):
                admissibleperm = self.admissible_perms[i]
                counter = 0 # number of copies of densperm in admissibleperm
                for c in combs:
                    admissible_subperm = normalize([admissibleperm.perm[x] for x in c], densperm.N)
                    admissible_subflag = PermFlag(admissible_subperm)
                    if admissible_subflag == densperm:
                       counter += 1
                       
                densities[i] += coeff*Integer(counter)/numcombs

        return densities
              


    def add_assumption(self, lincomb, bound):
        """Assumption of the form [(perm, coef), ..., (perm, coef)] \geq bound."""

        # switch on 'assumptions mode'
        if self.assumptions is False:
            self.assumptions = True

        assumption_densities = [-Rational(bound) for x in range(len(self.admissible_perms))]

        for densperm, coeff in lincomb:
            densperm = PermFlag(densperm)
            if self.N < densperm.N: # if too big to fit, then skip
                continue
            
            combs = Combinations(range(self.N), densperm.N)
            numcombs = combs.cardinality()
            
            # for each admissible permutation, compute the density of densperm in it
            for i in range(len(self.admissible_perms)):
                admissibleperm = self.admissible_perms[i]
                counter = 0 # number of copies of densperm in admissible perm
                for c in combs:
                    admissible_subperm = normalize([admissibleperm.perm[x] for x in c], densperm.N)
                    admissible_subflag = PermFlag(admissible_subperm)
                    if admissible_subflag == densperm:
                        counter += 1
                       
                assumption_densities[i] += coeff*Integer(counter)/numcombs
                #sys.stdout.write("assumption_densities[i] = %s\n" % str(assumption_densities[i]))

        self.assumption_densities.append(assumption_densities)
            

    def _export_flag_products(self, file_location="FP.txt"):

        """Store flag products in a file for future use."""

        """try:"""
        f = open(file_location, 'w')
        
        for ti in range(len(self.types)):
            f.write("%d\n" % ti)
            for fp in self.flag_products[ti]:
                f.write("%d %d %d %d %d\n" % (fp[0], fp[1], fp[2], fp[3], fp[4]))
                
        f.close()
        """except:
            raise IOError("Something went wrong with opening your file.\n")
        """

    
    def _load_flag_products(self, file_location="FP.txt"):

        """Retrieve flag products from a file in store."""

        """try:"""
        f = open(file_location, 'r')

        tp = 0
        self.flag_products = [[] for x in self.types]
        for line in f:
            ln = line.strip()
            if len(ln) < 6: #then it's type index
                tp = int(ln)
            else:
                fp = [Integer(x) for x in ln.split(' ')]
                self.flag_products[tp].append(fp)
        """except:
            raise IOError("Something went wrong when reading the file.\n")
        """
        
        
    def _compute_flag_products(self):
        
        """Parallel implementation of flag products computations."""

        # construct filename for this problem
        
        forb_str = "_".join([fp._to_str() for fp in self.forbidden_perms])

        filename = "FP-N" + str(self.N)
        if len(forb_str)>0:
            filename +="-forb_"+forb_str+".txt"
        else:
            filename +=".txt"

        # check if we have had this problem before (or any one with same flag_products)
        full_path = os.path.join(os.getcwd(),"..","store", filename)
        print full_path
        
        if os.path.exists(full_path):

            sys.stdout.write("Loading flag_products from file: ")
            sys.stdout.flush()
            self._load_flag_products(full_path)
            sys.stdout.write("OK.\n")

        else:
            
            sys.stdout.write("Computing flag products.\n")
            sys.stdout.write("Flag products on these types are done:\n")
            sys.stdout.flush()

            self.flag_products = [[] for tp in self.types]
            
            pool = multiprocessing.Pool(NUMCPU)
            num_tasks = len(self.types)
            result = pool.imap_unordered(self._compute_t_flag_products, range(num_tasks))

            """
            # PROGRESS REPORTING
            while True:
                completed = result._index
                sys.stderr.write("\rProgress: %d out of %d." % (completed, num_tasks))
                if completed == num_tasks:
                    break
            """

            for ti,tprods in result:
                self.flag_products[ti] = tprods
            pool.close()
            pool.join()

            sys.stdout.write("all done.\n")
            sys.stdout.flush()
            
            sys.stdout.write("Writing flag_products into file: ")
            sys.stdout.flush()
            self._export_flag_products(full_path) # store for future use
            sys.stdout.write("OK.\n")
            
        
        
    def _compute_t_flag_products(self, ti):
        """For each type t, generate all quintuples (s, p1, p2, a, b):

        s: admissible permutation (its index)
        p1: flag1 index in self.flags[t]
        p2: flag2 index in self.flags[t]
        a/b: density of p1*p2 in s

        10 Dec 15: for N=6, runs 4min28s
        """
        
        Nset = range(self.N)
        
        num_types = len(self.types)
        num_admissible_perms = len(self.admissible_perms)
        range_num_admissible_perms = range(num_admissible_perms)

        t = self.types[ti]
        m = int((self.N + t.N)/2) # flag order
        t_products = list() # will hold products of flags on type t
        # p(F1,F2,A) = 1/|\Theta|*\sum_{\theta} p(F1,F2,theta,A)
        denom = binomial(self.N,t.N)*binomial(self.N-t.N,m-t.N) 
        
        bigTheta = Combinations(self.N,t.N)#.list()
        
        num_ti_flags = len(self.flags[ti])
        range_num_ti_flags = range(num_ti_flags)
        
        for permi in range_num_admissible_perms:
                
            perm = self.admissible_perms[permi]
            
            # count number of times when F1xF2 == (perm,theta)  // over all theta
            counter =  [[0 for fi in range_num_ti_flags] for fj in range_num_ti_flags]
            
            for S0 in bigTheta: # S0 := im(theta)
                S0vals_in_perm = [perm.perm[x] for x in S0]
                nonimtheta = [x for x in Nset if x not in S0]
                S1_possibilities = Combinations(nonimtheta, m-t.N)
                    
                # condition 1
                if not PermFlag(normalize(S0vals_in_perm, t.N)) == t:
                    continue # go to next theta
                    
                # given cond1, how many times does (cond2 && cond3) hold for each fi,fj from ti-flags
                
                for S1 in S1_possibilities:
                    S1S0 = S1+S0
                    S1S0.sort()
                    subperm1 = [perm.perm[x] for x in S1S0]
                    subperm1_t = [subperm1.index(perm.perm[x]) for x in S0]
                    subperm2 = [perm.perm[x] for x in [x for x in Nset if x not in S1]]
                    subperm2_t = [subperm2.index(perm.perm[x]) for x in S0]
                    
                    perm1flag = PermFlag(normalize(subperm1, m), subperm1_t)
                    
                    toperm = normalize(subperm2,m)
                    perm2flag = PermFlag(toperm, subperm2_t)
                    
                    for fi in range_num_ti_flags:
                        for fj in range(fi, num_ti_flags): # only need to do distinct flags
                            
                            flg1 = self.flags[ti][fi]
                            flg2 = self.flags[ti][fj]
                            
                            # condition 2
                            if not perm1flag == flg1:
                                continue # go to the next S1
                                
                            # condition 3
                            if perm2flag == flg2:
                                counter[fi][fj] += 1
                        # end fj
                    # end fi
                # end S1
                
            # end S0 (i.e. theta)
            for fi in range_num_ti_flags:
                for fj in range_num_ti_flags:
                    if counter[fi][fj] > 0: # save to t_products
                        t_products.append([permi, fi, fj, counter[fi][fj], denom])

        sys.stdout.write("%d, " % ti)
        sys.stdout.flush()
        
        return ti,t_products

        

    def _write_sdp_file(self, filename="sdp.dat-s"):

        """Write SDP input file.

        INPUT:
        - filename: name of the file, extension should be .dat-s
        """

        
        sys.stdout.write("Writing SDP input file: ")
        sys.stdout.flush()

        self._sdp_file = filename
 
        ass_on = 0
        if self.assumptions:
            ass_on = 1
        
        num_types = len(self.types)
        num_constraints = len(self.admissible_perms) + 1 # need a constant too, need assumptions c_1+...+c_k = 1 constraint
        num_assumptions = len(self.assumption_densities)
        RHSvector = [0.0 for x in range(num_constraints)]
        RHSvector[-1] = 1.0 # const == 1 


        num_blocks = 1 + num_types + 2 + ass_on # (delta block) and (slacks block) and (\nu-constant block) and (1 if self.assumptions == True, 0 otherwise)
        block_sizes = [1]+[len(self.flags[i]) for i in range(num_types)]+[-len(self.admissible_perms), -1]
        if self.assumptions:
            block_sizes += [-num_assumptions]
        block_sizes_str = ' '.join([str(x) for x in block_sizes])
        
        sdpfile = open(filename, 'w')

        sdpfile.write(str(num_constraints)+"\n")
        sdpfile.write(str(num_blocks)+"\n")
        sdpfile.write(' '.join([str(x) for x in block_sizes])+"\n")
        sdpfile.write(' '.join([str(x) for x in RHSvector])+"\n")

        # write slacks and -delta to the LHS
        for i in range(len(self.admissible_perms)):
            sdpfile.write(str(i)+" 1 1 1 -1.0\n")
            sdpfile.write(str(i+1)+" "+str(1+num_types+1)+" "+str(i+1)+" "+str(i+1)+" 1.0\n")
        sdpfile.write(str(len(self.admissible_perms))+" 1 1 1 -1.0\n")

        # set constant to 1
        sdpfile.write(str(num_constraints)+" "+str(1+num_types+1+1)+" 1 1 1.0\n")
        
        # write densities to the LHS        
        for i in range(len(self.admissible_perms)):
            dens = self.densities[i]
            if dens != 0: # i is perm index
                sdpfile.write(str(i+1)+" "+str(num_blocks-ass_on)+" 1 1 "+str(dens.numerical_approx(digits=64))+"\n")

        # write assumptions
        if self.assumptions:
            # write them down
            for j in range(num_assumptions):
                for i in range(len(self.admissible_perms)):
                    ass_densi = self.assumption_densities[j][i]
                    if ass_densi != 0: # i is perm index
                        sdpfile.write(str(i+1)+" "+str(num_blocks)+" "+str(j+1)+" "+str(j+1)+" "+str(ass_densi.numerical_approx(digits=64))+"\n")
            """
            # write c_1+...+c_k = 1 constraint
            for j in range(num_assumptions):
                sdpfile.write(str(num_constraints)+" "+str(num_blocks)+" "+str(j+1)+" "+str(j+1)+" 1.0\n")
            """ 

        # write flag products to the LHS
        for i in range(num_types): # i: type
            for prod_entry in self.flag_products[i]:
                perm_index = prod_entry[0]+1 # 0: perm
                f1_index = prod_entry[1]+1   # 1: F1
                f2_index = prod_entry[2]+1   # 2: F2
                val = prod_entry[3]/prod_entry[4] # value
                
                sdpfile.write(str(perm_index)
                              +" "+str(1+i+1)
                              +" "+str(f1_index)
                              +" "+str(f2_index)
                              +" "+str(val.numerical_approx(digits=64))+"\n")

        sdpfile.close()
        
        sys.stdout.write("OK.\n")
        sys.stdout.flush()

        return

    def solve_sdp(self, outfile="sdp.out"):

        """ Ignore input for now, later can use it."""
        
        self._write_sdp_file()

        os.system("csdp sdp.dat-s sdp.out")
        

    def analyze_sdp_output(self, solfile="sdp.out"):

        """Take data from the output file created by the SDP solver."""
        
        sf = open(solfile, 'r')

        self._sdp_Q_matrices = [[[0 for x in self.flags[i]] for y in self.flags[i]] for i in range(len(self.types))]
        self._sdp_slacks = [0 for x in self.admissible_perms]
        
        for line in sf:
            l = line.strip().split(' ')

            # only process dual 
            if l[0] == '2':
                
                # clean line
                q = int(l[1]) # index of Q (+2)
                x = int(l[2])-1 # x coordinate in Q
                y = int(l[3])-1 # y coord in Q
                v = float(l[4]) # value at Q[x,y]
                
                # process line
                if q == 1:
                    self._sdp_bound = v

                if 1 < q < len(self.types)+2:
                    self._sdp_Q_matrices[q-2][x][y] = v
                    self._sdp_Q_matrices[q-2][y][x] = v

                if q == len(self.types)+2:
                    self._sdp_slacks[x] = v

        self._sdp_Q_matrices = [matrix(Q) for Q in self._sdp_Q_matrices]

        