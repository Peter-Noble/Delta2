import sympy
from sympy.physics.vector import *

x, y = sympy.symbols("x y")

Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Px, Py = sympy.symbols('ax ay bx by cx cy dx dy px py')
A, B, C, D, E, F, G, H, P, Pt = sympy.symbols('A B C D E F G H P Pt')
u, v = sympy.symbols("u v")

A = Ax * x + Ay * y
B = Bx * x + By * y
C = Cx * x + Cy * y
D = Dx * x + Dy * y
Pt = Px * x + Py * y

E = (1-v)*A + v*B
F = (1-v)*D + v*C
G = (1-u)*A + u*D
H = (1-u)*B + u*C

P = (1-v)*G + v*H
#P = (1-u)*E + u*F

print(P)

eq = P - Pt

print(eq)

eqx = eq.subs([(x, 1), (y, 0)])
eqy = eq.subs([(x, 0), (y, 1)])

print(eqx)
print(eqy)

result = sympy.solve([eqx, eqy], (u, v))

print(result)
