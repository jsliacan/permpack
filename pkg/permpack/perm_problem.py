import os, sys, time, multiprocessing, numpy, json

from numpy import array

from sage.misc.misc import SAGE_TMP
from sage.combinat.permutation import Permutation, Permutations, PatternAvoider
from sage.combinat.combination import Combinations
from sage.combinat.subset import Subsets, Set
from sage.arith.all import factorial, binomial
from sage.rings.all import Integer, Rational, QQ
from sage.doctest.util import Timer
from sage.misc.functional import numerical_approx
from sage.matrix.all import matrix, diagonal_matrix
from sage.modules.free_module_element import vector

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
        sys.stdout.write("Generating admissible permutations... ")
        sys.stdout.flush()
        self.admissible_perms = self._generate_admissible_perms()
        sys.stdout.write("\033[32mOK\033[m.\n")
        sys.stdout.flush()
        
        # TYPES
        sys.stdout.write("Generating types... ")
        sys.stdout.flush()
        self.types = self._generate_types()
        sys.stdout.write("\033[32mOK\033[m.\n")

        self._info = len(self.types)

        # FLAGS
        sys.stdout.write("Generating flags... ")
        sys.stdout.flush()
        # generate flags
        self.flags = self._generate_flags()
        sys.stdout.write("\033[32mOK\033[m.\n")

        # COMPUTING DENSITIES
        sys.stdout.write("Expressing density pattern as a linear combination of admissible permutations... ")
        sys.stdout.flush()
        self.densities = self._compute_densities()
        sys.stdout.write("\033[32mOK\033[m.\n")
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

        sys.stdout.write("Generating flags (SLOW VERSION).\n")
        sys.stdout.flush()
        
        flags = list()

        sys.stdout.write("Generating flags... ")
        
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

        sys.stdout.write("\033[32mOK\033[m.\n")
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

            sys.stdout.write("Loading flag_products from file... ")
            sys.stdout.flush()
            self._load_flag_products(full_path)
            sys.stdout.write("\033[32mOK\033[m.\n")

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

            sys.stdout.write("all done. \033[32mOK\033[m.\n")
            sys.stdout.flush()
            
            sys.stdout.write("Writing flag_products into file... ")
            sys.stdout.flush()
            self._export_flag_products(full_path) # store for future use
            sys.stdout.write("\033[32mOK\033[m.\n")

            
        
        
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

        
        sys.stdout.write("Writing SDP input file... ")
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
        
        sys.stdout.write("\033[32mOK\033[m.\n")
        sys.stdout.flush()

        return

    def solve_sdp(self, outfile="sdp.out", solver="csdp"):

        self._write_sdp_file()
        sys.stdout.write("Solving SDP problem...\n")
        sys.stdout.flush()

        if solver == "csdp":
            self._solver = "CSDP"
            os.system("csdp sdp.dat-s sdp.out")
        elif solver == "sdpa_dd":
            self._solver = "SDPA_DD"
            os.system("sdpa_dd -ds sdp.dat-s -o sdp.out")
        else:
            raise NotImplementedError

        sys.stdout.write("Finished. \033[32mOK\033[m.\n")

    def exactify(self, solfile="sdp.out", recognition_precision=10, rounding_precision=10e20):

        """
        Round matrices to have entries in a given field, so far only
        QQ. Does not accommodate assumptions.


        INPUT:

        - solfile: the name of the solution file, as a string

        - recognition_precision: number of decimal digits to round to 0

        - rounding_precison: denominator in rounding: round(a/rounding_precision)*rounding_precision

        EXAMPLE:

        sage: p = PermProblem(4, density_perm="132")
        sage: p.solve()
        sage: p.exactify(rounding_precision=10^10)
        sage: p.bound
        """

        sys.stdout.write("Transforming floating-point matrices to rational matrices...\n")
        sys.stdout.flush()
        
        num_types = len(self.types)
        num_flags = [len(self.flags[ti]) for ti in range(num_types)]
        
        
        # take data from the output file created by the SDP solver

        try:
            sf = open(solfile, 'r')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise

        self._sdp_Q_matrices = [[[0 for x in self.flags[i]] for y in self.flags[i]] for i in range(len(self.types))]
        self._sdp_slacks = [0 for x in self.admissible_perms]

        self._exact_Q_matrices = list()

        # solver: CSDP
        if self._solver == 'CSDP':
            sys.stdout.write("Reading output of the CSDP solver... ")
            sys.stdout.flush()
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

            sys.stdout.write("\033[32mOK\033[m.\n")
            
        # solver: SDPA-DD
        elif self._solver == 'SDPA_DD':

            sys.stdout.write("Reading output of the SDPA-DD solver... ")
            sys.stdout.flush()
            
            obj_val = None
            begin = False
            q = -1
            r = 0
            
            for line in sf:

                if line[:12] == "objValPrimal":
                    self._sdp_bound = float(line.strip().split(' ')[-1])
                    continue
                    
                # do nothing until yMat is found
                if line[:6] == "yMat =":
                    begin = True
                    continue

                # skip first/last bracket 
                if line.strip() == "{" or line.strip() == "}":
                    continue
                
                # skip the line block with obj value
                if begin and line[0] == "{" and obj_val == None:
                    
                    obj_val = self._sdp_bound
                    continue

                if begin and line[0] == "{":
                    
                    q += 1
                    r = 0
                    line = line.strip()
                    while line[0] in [' ','{', ',']:
                        line = line[1:]
                    while line[-1] in [' ','}', ',']:
                        line = line[:-1]
                    qrow = line.split(',')
                    
                    if q < len(self.types): # still reading in Q matrixes

                        for i in range(len(qrow)):
                            self._sdp_Q_matrices[q][r][i] = float(qrow[i])
                            
                    else: # then we are already reading in slacks

                        for i in range(len(qrow)):
                            self._sdp_slacks[i] = float(qrow[i])

                        break # we are finished.
                
                elif begin:
                    
                    r += 1
                    line = line.strip()
                    while line[0] in [' ', '{', ',']:
                        line = line[1:]
                    while line[-1] in [' ', '}', ',']:
                        line = line[:-1]

                    qrow = line.split(',')
                    for i in range(len(qrow)):
                        self._sdp_Q_matrices[q][r][i] = float(qrow[i])
                else:
                    continue
                
            sys.stdout.write("\033[32mOK\033[m.\n")
                        
        else:
            print "Your solver's output cannot be processed by Permpack. Sorry!"
            return

        
        self._sdp_Q_matrices = [matrix(Q) for Q in self._sdp_Q_matrices]

        
        # ROUNDING
        sys.stdout.write("Rounding Q matrices... ")
        sys.stdout.flush()
        
        for ti in range(num_types):
            Qi = self._sdp_Q_matrices[ti]

            # preprocess matrices (to avoid neg evalues)
            Lsdp = numpy.linalg.cholesky(Qi)


            L = matrix(QQ, num_flags[ti], num_flags[ti], sparse=True)
            for i in range(num_flags[ti]):
                for j in range(num_flags[ti]):
                    L[i,j] = round(Lsdp[i,j]*rounding_precision)/rounding_precision
            Qexact = L*L.transpose()
            self._exact_Q_matrices.append(Qexact)

        sys.stdout.write("\033[32mOK\033[m.\n")
        
        # FORM LINEAR SYSTEM OF EQUATIONS

        sys.stdout.write("Computing exact bound...\n")
        sys.stdout.flush()
            
        self._exact_bounds = [0 for x in self.admissible_perms]

        for ti in range(len(self.types)):
            for fp in self.flag_products[ti]:
                aij = Integer(fp[3])/Integer(fp[4])
                if fp[1] == fp[2]:
                    self._exact_bounds[fp[0]] += aij*self._exact_Q_matrices[ti][fp[1]][fp[2]]
                else:
                    self._exact_bounds[fp[0]] += aij*self._exact_Q_matrices[ti][fp[1]][fp[2]]*2

        for di in range(len(self.densities)):
            self._exact_bounds[di] += self.densities[di]
            
        self._exact_bound = max(self._exact_bounds)
        sys.stdout.write("\033[31m[WARNING] If you used assumptions, those are not included in the bounds! Matrices have still been rounded correctly. This will be added later.\033[m\n")
        sys.stdout.write("\033[32m[OK]   \033[mDone. Exact bound is roughly %.10f. Access it from self._exact_bound.\n"% self._exact_bound)
        sys.stdout.flush()
            
            

    """
        # transforming matrices

        self._Qdash_matrices = list()
        self._R_matrices = list()
        self._Rdash_matrices = list()
        self._Qexact_matrices = list()

        k = 0
        for Q in self._sdp_Q_matrices:

            sys.stdout.write("Transforming matrix Q for type %d\n" %k)
            sys.stdout.flush()

            
            sys.stdout.write("Finding eigen-stuff for matrix Q for type %d\n" %k)
            sys.stdout.flush()
            
            estuff = Q.eigenvectors_right()
            evectors = [x[1][0] for x in estuff] # right eigenvectors
            evals = [x[0] for x in estuff]

            sys.stdout.write("Computing zero-space for matrix Q for type %d\n" %k)
            sys.stdout.flush()
            
            zerospace_indices = list()
            zerospace_vectors = list()
            ones_indices = list()
            for i in range(len(evals)):
                if abs(evals[i]) < 1/recognition_precision: # will round these to 0
                    zerospace_indices.append(i)
                    evals[i] = 0
            
            # subtract vec from all other evectors so they are 0 where vec is 1
            for i in zerospace_indices:
                # make max vec_i entry equal 1
                maxvec = max(evectors[i])
                minvec = min(evectors[i])
                if abs(maxvec) > abs(minvec):
                    mindex = list(evectors[i]).index(maxvec)
                    mvec = maxvec
                else:
                    mindex = list(evectors[i]).index(minvec)
                    mvec = minvec
                veci = vector([x/mvec for x in evectors[i]])
                evectors[i] = veci
                ones_indices.append(mindex)

                # subtract vec_i from other vectors
                for j in zerospace_indices:
                    if j != i:
                        c = evectors[j][mindex]
                        veci_altered = vector([x*c for x in evectors[i]])
                        evectors[j] = evectors[j]-veci_altered

            for i in zerospace_indices:
                zerospace_vectors.append(evectors[i])
                
            # complete zerospace to basis for the whole space
            d = Q.dimensions()[0] # square matrix
            nonzerospace_basis = list()
            free_indices = range(d)
            for i in ones_indices:
                free_indices.remove(i)
                
            for i in range(d-len(zerospace_indices)):
                vec = [0 for x in range(d)]
                vec[free_indices[i]] = 1
                nonzerospace_basis.append(vec)

            sys.stdout.write("Computing R matrices for matrix Q for type %d\n" %k)
            sys.stdout.flush()
            
            # computing R matrices
            Rt = matrix(zerospace_vectors + nonzerospace_basis)
            Rtflat = Rt.list()
            for i in range(d*d):
                   Rtflat[i] = round(Rtflat[i]*rounding_precision)/rounding_precision # rational now
            R = matrix(QQ,d,d,Rtflat, sparse=True).transpose()
            self._R_matrices.append(R)

            sys.stdout.write("Computing Qdash matrices for Q for type %d\n" %k)
            sys.stdout.flush()
            
            # computing Qdash matrices
            Qdash = R.transpose()*Q*R
            # deleting rows/columns that will be multiplied by 0 eigenvectors in R
            Qdash = Qdash.delete_columns(range(len(zerospace_indices)))
            Qdash = Qdash.delete_rows(range(len(zerospace_indices)))
            qd = Qdash.dimensions()[0] # square mat
            Qdashflat = Qdash.list()
            for i in range(qd):
                Qdashflat[i] = round(Qdashflat[i]*rounding_precision)/rounding_precision # rational now
            Qdash = matrix(QQ,qd,qd,Qdashflat)
            #Qdash.delete_rows(range(len(zerospace_indices)))
            #Qdash.delete_columns(range(len(zerospace_indices)))
            self._Qdash_matrices.append(Qdash)

            sys.stdout.write("Computing Rdash matrices for Q for type %d\n" %k)
            sys.stdout.flush()
            start = time.time()
            # computing Rdash matrices
            # Rdash = R.inverse without the rows that correspond to 0 eigenvectors
            Rdash = R.inverse()
            end = time.time()
            sys.stdout.write("time duration: %f\n" %(end-start))
            Rdash = Rdash.delete_rows(range(len(zerospace_indices)))
            self._Rdash_matrices.append(Rdash)
            

            # Q = Rdash^T*Qdash*Rdash
            Qexact = Rdash.transpose()*Qdash*Rdash
            self._Qexact_matrices.append(Qexact)

            k += 1
        # matrices done

    """            

    def write_certificate(self, outfile="certificate.js"):

        """
        Write problem data to file so the bound can be verified. In particular,
        write problem description, admissible permutations, types, flags, L-matrices.

        """

        sys.stdout.write("Writing certificate into file... \n")
        sys.stdout.flush()

        density_str = ""
        i = 0

        if len(self.density_pattern) == 1 and self.density_pattern[0][1] == 1:
            density_str = str(self.density_pattern[0][0])
        else:
            for perm,coef in self.density_pattern:
                if i == 0:
                    density_str += str(coef)+"*"+str(perm)
                else:
                    density_str += " + "+str(coef)+"*"+str(perm)
                i += 1

        if self.forbidden_perms:
            problem_description = 'Maximize the density of '+ density_str+ " in "+ str(self.forbidden_perms)+"-free permutations."
        else:
            problem_description = 'Maximize the density of '+ density_str+"."


                
        data = {
            'problem': problem_description,
            'order of admissible permutations': int(self.N),
            'admissible permutations': [str(x) for x in self.admissible_perms],
            'types': [str(x) for x in self.types],
            'flags': [[str(x) for x in t_flags] for t_flags in self.flags],
            'bound': [str(x) for x in self._exact_bound]
        }
            
        
        with open(outfile, 'w') as certfile:
            json.dump(data, certfile, sort_keys=False, indent=4)

        sys.stdout.write("\033[32m[OK]   \033[mCertificate written successfully to \033[32m%s\033[m.\n" % outfile)
