from permpack.all import *

p = PermProblem(6, density_pattern="132")
p.solve_sdp()
p.analyze_sdp_output()
