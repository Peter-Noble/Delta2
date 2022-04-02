from sympy import *
from sympy.codegen.ast import Assignment, FunctionPrototype, FunctionDefinition, real, Return, Type, Declaration, Variable
from sympy.abc import n

gradient_type = Type("gradient")
float_type = Type("float")
vec_type = Type("vec3")

def dot2(x: Matrix):
    return x.dot(x)
