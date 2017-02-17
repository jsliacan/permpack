from permpack.all import *

p = PermProblem(6, density_pattern="1342")
p.solve_sdp(solver="csdp")
p.exactify(rounding_precision=10^20)

