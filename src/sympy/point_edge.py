from common import *

a, b, c, d, e, f, g, h, i, x, y, z = symbols("a b c d e f g h i x y z", real=True)
eps = symbols("eps", real=True)

A = Matrix([[a, b, c]])
B = Matrix([[d, e, f]])
C = Matrix([[g, h, i]])
X = Matrix([[x, y, z]])

AB = B - A
AC = C - A
AX = X - A

a00 = dot2(AB)
a01 = AB.dot(AC)
a11 = dot2(AC)
b0 = -AX.dot(AB)
b1 = -AX.dot(AC)

# E01
t0 = -b0/a00
t1 = 0

closest = A + t0 * AB + t1 * AC
distVec = X - closest
dist2 = dot2(distVec)
dist = sqrt(dist2)

adj = dist - eps

p = Piecewise((0, dist - eps < 0), (2*adj*adj, True))

vec_type = Type("Eigen::Vector3d")
vec_ref_type = Type("Eigen::Vector3d&")
void_type = Type("void")
float_type = Type("float")

params = [Variable("A", type=vec_type), Variable("B", type=vec_type), Variable("X", type=vec_type),
          Variable("eps", type=float_type),
          Variable("Ad", type=vec_ref_type), Variable("Bd", type=vec_ref_type)]
bdy = [
    "float a = A.x()",
    "float b = A.y()",
    "float c = A.z()",
    "float d = B.x()",
    "float e = B.y()",
    "float f = B.z()",
    "float x = X.x()",
    "float y = X.y()",
    "float z = X.z()",
]

for i, var in zip(["x()", "y()", "z()"], [a, b, c]):
    bdy.append("Ad.{} = {}".format(i, ccode(diff(p, var))))

for i, var in zip(["x()", "y()", "z()"], [d, e, f]):
    bdy.append("Bd.{} = {}".format(i, ccode(diff(p, var))))

pt_pt_grad = FunctionDefinition(return_type=void_type, name="pt_edge_gradient", parameters=params, body=bdy)
print(ccode(pt_pt_grad))
