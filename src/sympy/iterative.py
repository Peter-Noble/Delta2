import sympy
from sympy.printing.c import C11CodePrinter

a, b, c, d = sympy.symbols('a[t] b[t] c[t] d[t]')
a_update, b_update, c_update, d_update = sympy.symbols('a_update[t] b_update[t] c_update[t] d_update[t]')

alpha = sympy.symbols('alpha')
C_dontdividebyzero = sympy.symbols('C_dontdividebyzero')

# 1-to-n
"""xa11, xa12, xa13 = sympy.symbols('x_a11 x_a12 x_a13')  # pt 1 x, y, z
xa21, xa22, xa23 = sympy.symbols('x_a21 x_a22 x_a23')
xa31, xa32, xa33 = sympy.symbols('x_a31 x_a32 x_a33')"""

# 1-to-1-n-times
xa11, xa12, xa13 = sympy.symbols('x_a11[t] x_a12[t] x_a13[t]')  # pt 1 x, y, z
xa21, xa22, xa23 = sympy.symbols('x_a21[t] x_a22[t] x_a23[t]')
xa31, xa32, xa33 = sympy.symbols('x_a31[t] x_a32[t] x_a33[t]')

xb11, xb12, xb13 = sympy.symbols('x_b11[t] x_b12[t] x_b13[t]')
xb21, xb22, xb23 = sympy.symbols('x_b21[t] x_b22[t] x_b23[t]')
xb31, xb32, xb33 = sympy.symbols('x_b31[t] x_b32[t] x_b33[t]')

distance1 = (xa11 + a * (xa21 - xa11) + b * (xa31 - xa11)) - (xb11 + c * (xb21 - xb11) + d * (xb31 - xb11))
distance2 = (xa12 + a * (xa22 - xa12) + b * (xa32 - xa12)) - (xb12 + c * (xb22 - xb12) + d * (xb32 - xb12))
distance3 = (xa13 + a * (xa23 - xa13) + b * (xa33 - xa13)) - (xb13 + c * (xb23 - xb13) + d * (xb33 - xb13))

"""penalty = sympy.Max(0, -a) + sympy.Max(0, -b) + sympy.Max(0, -c) + sympy.Max(0, -d) \
               + sympy.Max(0, a - 1) + sympy.Max(0, b - 1) + sympy.Max(0, c - 1) \
               + sympy.Max(0, d - 1) # + sympy.Max(0, a + b - 1) + sympy.Max(0, c + d - 1)"""
penalty = sympy.Max(0, -a) + sympy.Max(0, -b) + sympy.Max(0, -c) + sympy.Max(0, -d) \
               + sympy.Max(0, a + b - 1) + sympy.Max(0, c + d - 1)

#J = 0.5 * (distance1 ** 2 + distance2 ** 2 + distance3 ** 2 + alpha * penalty ** 2)
J = 0.5 * (distance1 ** 2 + distance2 ** 2 + distance3 ** 2)

# J = 0.5 * (distance1 ** 2 + distance2 ** 2 + distance3 ** 2)

dJ_a = sympy.simplify(sympy.diff(J, a))
dJ_b = sympy.simplify(sympy.diff(J, b))
dJ_c = sympy.simplify(sympy.diff(J, c))
dJ_d = sympy.simplify(sympy.diff(J, d))

ddJ_aa = sympy.simplify(sympy.diff(dJ_a, a))
ddJ_bb = sympy.simplify(sympy.diff(dJ_b, b))
ddJ_cc = sympy.simplify(sympy.diff(dJ_c, c))
ddJ_dd = sympy.simplify(sympy.diff(dJ_d, d))

update_a = a - 1.0 / (C_dontdividebyzero + ddJ_aa) * dJ_a
update_b = b - 1.0 / (C_dontdividebyzero + ddJ_bb) * dJ_b
update_c = c - 1.0 / (C_dontdividebyzero + ddJ_cc) * dJ_c
update_d = d - 1.0 / (C_dontdividebyzero + ddJ_dd) * dJ_d

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


update_a = update_a.replace(sympy.DiracDelta(wild), 0)
update_b = update_b.replace(sympy.DiracDelta(wild), 0)
update_c = update_c.replace(sympy.DiracDelta(wild), 0)
update_d = update_d.replace(sympy.DiracDelta(wild), 0)

printer = CustomPrinter()

ccode_a = "a_new[t] = " + printer.doprint(update_a).splitlines()[-1] + ";"
ccode_b = "b_new[t] = " + printer.doprint(update_b).splitlines()[-1] + ";"
ccode_c = "c_new[t] = " + printer.doprint(update_c).splitlines()[-1] + ";"
ccode_d = "d_new[t] = " + printer.doprint(update_d).splitlines()[-1] + ";"

print("""
auto Heaviside = [](double x)->double {return x<0 ? 0.0 : 1.0;};
const int Iterations    = 4;
const int TrianglePairs = 8;
for (int i=0; i<Iterations; i++)
#pragma omp simd
for (int t=0; t<TrianglePairs; t++) {
""")
print(ccode_a.replace("1.0*", "1.0f*"))
print(ccode_b.replace("1.0*", "1.0f*"))
print(ccode_c.replace("1.0*", "1.0f*"))
print(ccode_d.replace("1.0*", "1.0f*"))
print("""
}
""")
