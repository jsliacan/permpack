from permpack.all import *

p = PermProblem(8, forbid=["312","231"], density_pattern="1324")
p.solve_sdp(show_output=True, solver="csdp")
p.exactify(rounding_precision=10^10, recognition_precision=6)
