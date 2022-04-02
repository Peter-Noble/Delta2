from common import *


a, b, c, x, y, z = symbols("a b c x y z", real=True)
eps = symbols("eps", real=True)

A = Matrix([[a, b, c]])
X = Matrix([[x, y, z]])

AX = X - A

dist = sqrt(dot2(AX))

p = Piecewise((0, dist - eps < 0), (2*(sqrt(dot2(AX)) - eps)**2, True))

vec_type = Type("Eigen::Vector3d")
vec_ref_type = Type("Eigen::Vector3d&")
void_type = Type("void")
float_type = Type("float")

params = [Variable("A", type=vec_type), Variable("X", type=vec_type), Variable("eps", type=float_type), Variable("Ad", type=vec_ref_type)]

bdy = [
    "float a = A.x()",
    "float b = A.y()",
    "float c = A.z()",
    "float x = X.x()",
    "float y = X.y()",
    "float z = X.z()",
]

for i, var in zip(["x()", "y()", "z()"], [a, b, c]):
    bdy.append("Ad.{} = {}".format(i, ccode(diff(p, var))))

pt_pt_grad = FunctionDefinition(return_type=void_type, name="pt_pt_gradient", parameters=params, body=bdy)
print(ccode(pt_pt_grad))
