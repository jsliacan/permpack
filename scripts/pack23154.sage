"""
Exact bound of ~0.16039352 on N=6.
"""

from permpack.all import *

p = PermProblem(6, density_pattern="23154")
p.solve_sdp(solver="csdp")
p.exactify(rounding_precision=10^10, recognition_precision=6)
p.write_certificate("cert23154.js")
