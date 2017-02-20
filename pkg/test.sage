from permpack.all import *


p = PermProblem(3, density_pattern="132")
p.solve_sdp()
p.exactify()
p.write_certificate()