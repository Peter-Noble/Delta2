from common import *

a, b, c, d, e, f, g, h, l, x, y, z = symbols("a b c d e f g h l x y z", real=True)
eps = symbols("eps", real=True)

A = Matrix([[a, b, c]])
B = Matrix([[d, e, f]])
C = Matrix([[g, h, l]])
X = Matrix([[x, y, z]])

AB = B - A
AC = C - A
AX = X - A

N = AB.cross(AC)
D = N.dot(A)

dist = (N.dot(X) - D) / sqrt(dot2(N))
dist2 = dist*dist

p = Piecewise((0, dist - eps < 0), (2*(sqrt(dist2) - eps)**2, True))

vec_type = Type("Eigen::Vector3d")
vec_ref_type = Type("Eigen::Vector3d&")
void_type = Type("void")

params = [Variable("A", type=vec_type), Variable("B", type=vec_type), Variable("C", type=vec_type),
          Variable("X", type=vec_type), Variable("eps", type=float_type),
          Variable("Ad", type=vec_ref_type), Variable("Bd", type=vec_ref_type), Variable("Cd", type=vec_ref_type)]
bdy = [
    "float a = A.x()",
    "float b = A.y()",
    "float c = A.z()",
    "float d = B.x()",
    "float e = B.y()",
    "float f = B.z()",
    "float g = C.x()",
    "float h = C.y()",
    "float l = C.z()", # i is already taken by sympy
    "float x = X.x()",
    "float y = X.y()",
    "float z = X.z()",
]

for i, var in zip(["x()", "y()", "z()"], [a, b, c]):
    bdy.append("Ad.{} = {}".format(i, ccode(diff(p, var))))

for i, var in zip(["x()", "y()", "z()"], [d, e, f]):
    bdy.append("Bd.{} = {}".format(i, ccode(diff(p, var))))

for i, var in zip(["x()", "y()", "z()"], [g, h, l]):
    bdy.append("Cd.{} = {}".format(i, ccode(diff(p, var))))

pt_pt_grad = FunctionDefinition(return_type=void_type, name="pt_face_gradient", parameters=params, body=bdy)
print(ccode(pt_pt_grad))
