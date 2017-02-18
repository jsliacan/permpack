from permpack.all import *

p = PermProblem(7, density_pattern="1342")
p.solve_sdp(solver="csdp")
p.exactify(rounding_precision=10^20)

