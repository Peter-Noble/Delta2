import sympy
from sympy.physics.vector import *

N = ReferenceFrame('N')
t = sympy.symbols("t")

t_0_x, t_0_y, t_0_z = sympy.symbols("t0x t0y t0z")
t_1_x, t_1_y, t_1_z = sympy.symbols("t1x t1y t1z")
t_2_x, t_2_y, t_2_z = sympy.symbols("t2x t2y t2z")

p_start_x, p_start_y, p_start_z = sympy.symbols("psx psy psz")
p_end_x, p_end_y, p_end_z = sympy.symbols("pex pey pez")
p = sympy.symbols("p")

p = (p_start_x * N.x + p_start_y * N.y + p_start_z * N.z) * (1 - t) + (p_end_x * N.x + p_end_y * N.y + p_end_z * N.z) * t

t_0, t_1, t_2 = sympy.symbols("ts0 ts1 ts2")

t_0 = N.x * t_0_x + N.y * t_0_y + N.z * t_0_z
t_1 = N.x * t_1_x + N.y * t_1_y + N.z * t_1_z
t_2 = N.x * t_2_x + N.y * t_2_y + N.z * t_2_z

n = (t_1 - t_0).cross(t_2 - t_0)

AP = sympy.symbols("AP")
f = sympy.symbols("f")

AP = p - t_0

f = n.dot(AP)

print(f)
print("===")
print(f.subs({p_start_x:0, p_end_x:0, p_start_y:0, p_end_y: 0, p_start_z:0, p_end_y:0}))
print("===")

result = sympy.solve(f, t)

print(result)

print("===")

print(sympy.cse(result))
