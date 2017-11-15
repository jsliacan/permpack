from permpack.all import *

p = PermProblem(7, density_pattern="1342")
p.solve_sdp(show_output=True, solver="csdp")
p.exactify(rounding_precision=10^10, recognition_precision=6)