from flagmatic.all import *

F = ["3:12","3:1223", "3:122331"]
S = "4:1213142434"
p = OrientedGraphProblem(8, forbid_induced=F, density=S)
p.solve_sdp(show_output=True, solver="csdp")
p.make_exact(10^20)