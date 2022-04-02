import sympy


a0startx, a0starty, a0startz = sympy.symbols('a_start.A.x() a_start.A.y() a_start.A.z()')
a1startx, a1starty, a1startz = sympy.symbols('a_start.B.x() a_start.B.y() a_start.B.z()')
a0endx, a0endy, a0endz = sympy.symbols('a_end.A.x() a_end.A.y() a_end.A.z()')
a1endx, a1endy, a1endz = sympy.symbols('a_end.B.x() a_end.B.y() a_end.B.z()')

b0startx, b0starty, b0startz = sympy.symbols('b_start.A.x() b_start.A.y() b_start.A.z()')
b1startx, b1starty, b1startz = sympy.symbols('b_start.B.x() b_start.B.y() b_start.B.z()')
b0endx, b0endy, b0endz = sympy.symbols('b_end.A.x() b_end.A.y() b_end.A.z()')
b1endx, b1endy, b1endz = sympy.symbols('b_end.B.x() b_end.B.y() b_end.B.z()')

t, pa, pb = sympy.symbols("t pa pb")

# t0x, t0y, t0z = sympy.symbols('t0x[id] t0y[id] t0z[id]')
# t1x, t1y, t1z = sympy.symbols('t1x[id] t1y[id] t1z[id]')

tax = (a0startx * (1 - t) + a0endx * t) * (1 - pa) + (a1startx * (1 - t) + a1endx) * pa
tay = (a0starty * (1 - t) + a0endy * t) * (1 - pa) + (a1starty * (1 - t) + a1endy) * pa
taz = (a0startz * (1 - t) + a0endz * t) * (1 - pa) + (a1startz * (1 - t) + a1endz) * pa

tbx = (b0startx * (1 - t) + b0endx * t) * (1 - pb) + (b1startx * (1 - t) + b1endx) * pb
tby = (b0starty * (1 - t) + b0endy * t) * (1 - pb) + (b1starty * (1 - t) + b1endy) * pb
tbz = (b0startz * (1 - t) + b0endz * t) * (1 - pb) + (b1startz * (1 - t) + b1endz) * pb

dist = sympy.symbols('dist')

dist = sympy.sqrt((tax - tbx)**2 + (tay - tby)**2 + (taz - tbz)**2)

diff = sympy.diff(dist, t)
diff2 = sympy.diff(diff, t)
print(diff)
print(diff2)
