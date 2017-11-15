"""
Exact bound of ~0.0673094718 on N=6
"""

from permpack.all import *

p = PermProblem(6, density_pattern="231564")
p.solve_sdp(solver="csdp")
p.exactify(rounding_precision=10^10, recognition_precision=6)
p.write_certificate("cert231564.js")
