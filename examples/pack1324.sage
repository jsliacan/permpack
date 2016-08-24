from permpack.all import *

# we know the maximizer is layered -> can forbid things
p = PermProblem(7, density_pattern="1324", forbid=["231","312"])
p.solve_sdp()
