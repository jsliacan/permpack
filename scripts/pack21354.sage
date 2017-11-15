"""
Exact bound of ~0.165150672 on N=7

* can forbid 231 and 312 because we know that there's a layered maximiser
"""

from permpack.all import *

p = PermProblem(7, forbid=["231","312"], density_pattern="21354")
p.solve_sdp(solver="csdp")
p.exactify(rounding_precision=10^10, recognition_precision=6)
p.write_certificate("cert21354.js")