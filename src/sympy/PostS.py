import sympy
from sympy.physics.vector import *
from sympy import sqrt, log, Max, Min
 
# d0, d1, v0, v1, Meff, dt = sympy.symbols("d_0 d_1 v_0 v_1 M t")  # Distance along normal at t=0, v at t=0, current estimate of v at t=0, effective mass along contact normal, timestep size
# inner, outer = sympy.symbols("i o")
# s, vs, ds, delta_s, F_ext= sympy.symbols("s v_s d_s delta_s F_ext")  # The thing we're trying to solve for
# s_last, T, dc, x, F_last, F_total = sympy.symbols("s_l T d_c x F_last F_total")

d0, d1, v0, v1, Meff, dt = sympy.symbols("d0 d1 v0 v1 Meff dt")  # Distance along normal at t=0, v at t=0, current estimate of v at t=0, effective mass along contact normal, timestep size
inner, outer = sympy.symbols("inner outer")
s, vs, ds, delta_s, F_ext= sympy.symbols("s vs ds delta_s F_ext")  # The thing we're trying to solve for
d_start = sympy.symbols("d_start")
s_last, dc, x, F_last, F_total = sympy.symbols("s_last dc x F_last F_total")


def P(dist):
   return (inner - outer) / (dist - inner) - dist / (inner - outer)

def Fi(dist):
   return Max(0.0, (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner)))

def F(dist):
   return Fi(Max(Min(dist, outer), 0.0))

t = (1-s)*dt
final_pos = d1 + dc
move = ds + vs * t + 0.5 * (((F(ds) + F(d1 + dc)) / 2.0 + F_total) / Meff) * t * t

eq = final_pos - move

print("F")
print(str(F(x)).replace("**", "^").replace("Max", "max").replace("Min", "min"))

print("eq")
print(str(eq).replace("**", "^").replace("Max", "max").replace("Min", "min"))

print("grad")
print(sympy.diff(eq, dc))

# print("dc solve")
# print(str(sympy.solve(eq, dc)).replace("**", "^").replace("Max", "max").replace("Min", "min"))
