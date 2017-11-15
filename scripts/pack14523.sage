"""
Exact bound of ~0.1536498380 on N=6
"""

from permpack.all import *

p = PermProblem(6, density_pattern="14523")
p.solve_sdp(solver="csdp")
p.exactify(rounding_precision=10^10, recognition_precision=6)
p.write_certificate("cert14523.js")