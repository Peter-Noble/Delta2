import sympy

# https://stackoverflow.com/questions/35819142/calculate-a-2d-homogeneous-perspective-transformation-matrix-from-4-points-in-ma

Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Px, Py = sympy.symbols('ax ay bx by cx cy dx dy px py')

M = sympy.Matrix([
    [Ax, Ay, 1, 0, 0, 0, 0, 0],
    [Bx, By, 1, 0, 0, 0, 0, 0],
    [Cx, Cy, 1, 0, 0, 0, -Cx, -Cy],
    [Dx, Dy, 1, 0, 0, 0, -Dx, -Dy],
    [0, 0, 0, Ax, Ay, 1, 0, 0],
    [0, 0, 0, Bx, Cy, 1, -Bx, -By],
    [0, 0, 0, Cx, Cy, 1, -Cx, -Cy],
    [0, 0, 0, Dx, Dy, 1, 0, 0]
])

b = sympy.Matrix(8, 1, [0, 0, 1, 1, 0, 1, 1, 0])

c00, c01, c02, c10, c11, c12, c20, c21 = M.LUsolve(b)

T = sympy.Matrix([
    [c00, c01, c02],
    [c10, c11, c12],
    [c20, c21, 1]
])

print(T)

print(T**-1)
