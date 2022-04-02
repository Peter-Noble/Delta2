from common import *

vec_type = Type("Eigen::Vector3d")
vec_ref_type = Type("Eigen::Vector3d&")
void_type = Type("void")
float_type = Type("float")

a, b, c, d, e, f, g, h, l = symbols("a b c d e f g h l", real=True)

A = Matrix([[a, b, c]])
B = Matrix([[d, e, f]])
C = Matrix([[g, h, l]])

area = 0.5 * dot2((A-C).cross(B-C))

params = [Variable("A", type=vec_type), Variable("B", type=vec_type), Variable("C", type=vec_type),
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
    "float l = C.z()"
]

for i, var in zip(["x()", "y()", "z()"], [a, b, c]):
    bdy.append("Ad.{} = {}".format(i, ccode(diff(2 * area * area, var))))

for i, var in zip(["x()", "y()", "z()"], [d, e, f]):
    bdy.append("Bd.{} = {}".format(i, ccode(diff(2 * area * area, var))))

for i, var in zip(["x()", "y()", "z()"], [g, h, l]):
    bdy.append("Cd.{} = {}".format(i, ccode(diff(2 * area * area, var))))

area_grad = FunctionDefinition(return_type=void_type, name="area_gradient", parameters=params, body=bdy)
print(ccode(area_grad))

AB = B - A
BC = B - C
CA = C - A

degen = (1/(1/(dot2(AB)**n) + 1/(dot2(BC)**n) + 1/(dot2(CA)**n))) / (dot2(AB)**n + dot2(BC)**n + dot2(CA)**n)

bdy = [
    "const int n = 4",
    "float a = A[0]",
    "float b = A[1]",
    "float c = A[2]",
    "float d = B[0]",
    "float e = B[1]",
    "float f = B[2]",
    "float g = C[0]",
    "float h = C[1]",
    "float l = C[2]"
]

for i, var in zip(["x()", "y()", "z()"], [a, b, c]):
    bdy.append("Ad.{} = {}".format(i, ccode(diff(2 * degen * degen, var))))

for i, var in zip(["x()", "y()", "z()"], [d, e, f]):
    bdy.append("Bd.{} = {}".format(i, ccode(diff(2 * degen * degen, var))))

for i, var in zip(["x()", "y()", "z()"], [g, h, l]):
    bdy.append("Cd.{} = {}".format(i, ccode(diff(2 * degen * degen, var))))

degen_grad = FunctionDefinition(return_type=void_type, name="degen_gradient", parameters=params, body=bdy)
print(ccode(degen_grad))
