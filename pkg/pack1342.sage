from permpack.all import *

p = PermProblem(7, density_pattern="1342")
p.solve_sdp()
p.exactify(rounding_precision=10^20)
p.write_certificate("cert1342.js")