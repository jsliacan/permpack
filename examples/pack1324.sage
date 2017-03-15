from flagmatic.all import *

p = OrientedGraphProblem(8, forbid_induced=["3:12", "3:1223", "3:122331"], density="4:1213142434")
p.solve_sdp(solver="csdp")
p.make_exact(denominator=10^20)
p.write_certificate("cert1324.js")
