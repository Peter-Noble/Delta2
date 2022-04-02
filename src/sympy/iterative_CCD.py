import sympy
from sympy.printing.c import C11CodePrinter

a, b, c, d = sympy.symbols('a[id] b[id] c[id] d[id]')
a_update, b_update, c_update, d_update, t_update = sympy.symbols('a_update[id] b_update[id] c_update[id] d_update[id] t_update[id]')

alpha = sympy.symbols('alpha')
C_dontdividebyzero = sympy.symbols('C_dontdividebyzero')

t = sympy.symbols("t[id]")

# 1-to-n
"""xa11, xa12, xa13 = sympy.symbols('x_a11 x_a12 x_a13')  # pt 1 x, y, z
xa21, xa22, xa23 = sympy.symbols('x_a21 x_a22 x_a23')
xa31, xa32, xa33 = sympy.symbols('x_a31 x_a32 x_a33')"""

# 1-to-1-n-times
xa11_t0, xa12_t0, xa13_t0 = sympy.symbols('x_a11_t0[id] x_a12_t0[id] x_a13_t0[id]')  # pt 1 x, y, z
xa21_t0, xa22_t0, xa23_t0 = sympy.symbols('x_a21_t0[id] x_a22_t0[id] x_a23_t0[id]')
xa31_t0, xa32_t0, xa33_t0 = sympy.symbols('x_a31_t0[id] x_a32_t0[id] x_a33_t0[id]')

xa11_t1, xa12_t1, xa13_t1 = sympy.symbols('x_a11_t1[id] x_a12_t1[id] x_a13_t1[id]')
xa21_t1, xa22_t1, xa23_t1 = sympy.symbols('x_a21_t1[id] x_a22_t1[id] x_a23_t1[id]')
xa31_t1, xa32_t1, xa33_t1 = sympy.symbols('x_a31_t1[id] x_a32_t1[id] x_a33_t1[id]')

xb11_t0, xb12_t0, xb13_t0 = sympy.symbols('x_b11_t0[id] x_b12_t0[id] x_b13_t0[id]')
xb21_t0, xb22_t0, xb23_t0 = sympy.symbols('x_b21_t0[id] x_b22_t0[id] x_b23_t0[id]')
xb31_t0, xb32_t0, xb33_t0 = sympy.symbols('x_b31_t0[id] x_b32_t0[id] x_b33_t0[id]')

xb11_t1, xb12_t1, xb13_t1 = sympy.symbols('x_b11_t1[id] x_b12_t1[id] x_b13_t1[id]')
xb21_t1, xb22_t1, xb23_t1 = sympy.symbols('x_b21_t1[id] x_b22_t1[id] x_b23_t1[id]')
xb31_t1, xb32_t1, xb33_t1 = sympy.symbols('x_b31_t1[id] x_b32_t1[id] x_b33_t1[id]')

xa11, xa12, xa13 = sympy.symbols('x_a11[t] x_a12[t] x_a13[t]')  # itermediate vars
xa21, xa22, xa23 = sympy.symbols('x_a21[t] x_a22[t] x_a23[t]')
xa31, xa32, xa33 = sympy.symbols('x_a31[t] x_a32[t] x_a33[t]')

xb11, xb12, xb13 = sympy.symbols('x_b11[t] x_b12[t] x_b13[t]')
xb21, xb22, xb23 = sympy.symbols('x_b21[t] x_b22[t] x_b23[t]')
xb31, xb32, xb33 = sympy.symbols('x_b31[t] x_b32[t] x_b33[t]')

xa11 = xa11_t0 * (1 - t) + xa11_t1 * t
xa12 = xa12_t0 * (1 - t) + xa12_t1 * t
xa13 = xa13_t0 * (1 - t) + xa13_t1 * t
xa21 = xa21_t0 * (1 - t) + xa21_t1 * t
xa22 = xa22_t0 * (1 - t) + xa22_t1 * t
xa23 = xa23_t0 * (1 - t) + xa23_t1 * t
xa31 = xa31_t0 * (1 - t) + xa31_t1 * t
xa32 = xa32_t0 * (1 - t) + xa32_t1 * t
xa33 = xa33_t0 * (1 - t) + xa33_t1 * t

xb11 = xb11_t0 * (1 - t) + xb11_t1 * t
xb12 = xb12_t0 * (1 - t) + xb12_t1 * t
xb13 = xb13_t0 * (1 - t) + xb13_t1 * t
xb21 = xb21_t0 * (1 - t) + xb21_t1 * t
xb22 = xb22_t0 * (1 - t) + xb22_t1 * t
xb23 = xb23_t0 * (1 - t) + xb23_t1 * t
xb31 = xb31_t0 * (1 - t) + xb31_t1 * t
xb32 = xb32_t0 * (1 - t) + xb32_t1 * t
xb33 = xb33_t0 * (1 - t) + xb33_t1 * t

distance1 = (xa11 + a * (xa21 - xa11) + b * (xa31 - xa11)) - (xb11 + c * (xb21 - xb11) + d * (xb31 - xb11))
distance2 = (xa12 + a * (xa22 - xa12) + b * (xa32 - xa12)) - (xb12 + c * (xb22 - xb12) + d * (xb32 - xb12))
distance3 = (xa13 + a * (xa23 - xa13) + b * (xa33 - xa13)) - (xb13 + c * (xb23 - xb13) + d * (xb33 - xb13))

"""penalty = sympy.Max(0, -a) + sympy.Max(0, -b) + sympy.Max(0, -c) + sympy.Max(0, -d) \
               + sympy.Max(0, a - 1) + sympy.Max(0, b - 1) + sympy.Max(0, c - 1) \
               + sympy.Max(0, d - 1) # + sympy.Max(0, a + b - 1) + sympy.Max(0, c + d - 1)"""
penalty = sympy.Max(0, -a) + sympy.Max(0, -b) + sympy.Max(0, -c) + sympy.Max(0, -d) \
               + sympy.Max(0, a + b - 1) + sympy.Max(0, c + d - 1) \
                   + sympy.Max(0, -t) + sympy.Max(0, t - 1)

#J = 0.5 * (distance1 ** 2 + distance2 ** 2 + distance3 ** 2 + alpha * penalty ** 2)
J = 0.5 * (distance1 ** 2 + distance2 ** 2 + distance3 ** 2)

#J = penalty

# J = 0.5 * (distance1 ** 2 + distance2 ** 2 + distance3 ** 2)


#dJ_a = sympy.simplify(sympy.diff(J, a))
#dJ_b = sympy.simplify(sympy.diff(J, b))
#dJ_c = sympy.simplify(sympy.diff(J, c))
#dJ_d = sympy.simplify(sympy.diff(J, d))
dJ_t = sympy.simplify(sympy.diff(J, t))

print("=== dJ_t ====")
print(dJ_t)


#ddJ_aa = sympy.simplify(sympy.diff(dJ_a, a))
#ddJ_bb = sympy.simplify(sympy.diff(dJ_b, b))
#ddJ_cc = sympy.simplify(sympy.diff(dJ_c, c))
#ddJ_dd = sympy.simplify(sympy.diff(dJ_d, d))
ddJ_tt = sympy.simplify(sympy.diff(dJ_t, t))

print("=== ddJ_tt ====")
print(ddJ_tt)


#update_a = a - sympy.simplify(1.0 / (C_dontdividebyzero + ddJ_aa) * dJ_a)
#update_b = b - sympy.simplify(1.0 / (C_dontdividebyzero + ddJ_bb) * dJ_b)
#update_c = c - sympy.simplify(1.0 / (C_dontdividebyzero + ddJ_cc) * dJ_c)
#update_d = d - sympy.simplify(1.0 / (C_dontdividebyzero + ddJ_dd) * dJ_d)
#update_t = t - sympy.simplify(1.0 / (C_dontdividebyzero + ddJ_tt) * dJ_t)
update_t = t - 1.0 / (C_dontdividebyzero + ddJ_tt) * dJ_t

print("=== update_t ====")
print(update_t)

# update_a = a - alpha * dJ_a
# update_b = b - alpha * dJ_b
# update_c = c - alpha * dJ_c
# update_d = d - alpha * dJ_d
# update_t = t - alpha * dJ_t

# print("One update rule:")
# print("================")
# print("a: ")
# print(str(update_a))
# print(sympy.latex(update_a))
# print("b: ")
# print(str(update_b))
# print(sympy.latex(update_b))
# print("c: ")
# print(str(update_c))
# print(sympy.latex(update_c))
# print("d: ")
# print(str(update_d))
# print(sympy.latex(update_d))

# print("C code for one update (raw output):")
# print("======================")
wild = sympy.Wild("wild")


class CustomPrinter(C11CodePrinter):
    def _print_Pow(self, expr):
        if expr.exp == 2:
            return self._print("({0}) * ({0})".format(expr.base))
        else:
            return self._print(expr)


#update_a = update_a.replace(sympy.DiracDelta(wild), 0)
#update_b = update_b.replace(sympy.DiracDelta(wild), 0)
#update_c = update_c.replace(sympy.DiracDelta(wild), 0)
#update_d = update_d.replace(sympy.DiracDelta(wild), 0)
update_t = update_t.replace(sympy.DiracDelta(wild), 0)

printer = CustomPrinter()

#ccode_a = "a_new[id] = " + printer.doprint(update_a).splitlines()[-1] + ";"
#ccode_b = "b_new[id] = " + printer.doprint(update_b).splitlines()[-1] + ";"
#ccode_c = "c_new[id] = " + printer.doprint(update_c).splitlines()[-1] + ";"
#ccode_d = "d_new[id] = " + printer.doprint(update_d).splitlines()[-1] + ";"
ccode_t = "t_new[id] = " + printer.doprint(update_t).splitlines()[-1] + ";"

print("""
auto Heaviside = [](double x)->double {return x<0 ? 0.0 : 1.0;};
const int Iterations    = 4;
const int TrianglePairs = 8;
for (int i=0; i<Iterations; i++)
#pragma omp simd
for (int t=0; t<TrianglePairs; t++) {
""")
#print(ccode_a.replace("1.0*", "1.0f*"))
#print(ccode_b.replace("1.0*", "1.0f*"))
#print(ccode_c.replace("1.0*", "1.0f*"))
#print(ccode_d.replace("1.0*", "1.0f*"))
print(ccode_t.replace("1.0*", "1.0f*"))
print("""
}
""")
