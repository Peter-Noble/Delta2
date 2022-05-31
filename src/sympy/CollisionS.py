import sympy
from sympy.physics.vector import *
from sympy import sqrt, log
 
d0, v0, v1, Meff, dt = sympy.symbols("d_0 v_0 v_1 M t")  # Distance along normal at t=0, v at t=0, current estimate of v at t=0, effective mass along contact normal, timestep size
inner, outer = sympy.symbols("i o")
s, vs, ds, delta_s, F_ext= sympy.symbols("s v_s d_s delta_s F_ext")  # The thing we're trying to solve for
s_last = sympy.symbols("s_l")

# d0, v0, v1, Meff, dt = sympy.symbols("d0 v0 v1 Meff dt")  # Distance along normal at t=0, v at t=0, current estimate of v at t=0, effective mass along contact normal, timestep size
# inner, outer = sympy.symbols("inner outer")
# s, vs, ds, delta_s, F_ext= sympy.symbols("s vs ds delta_s F_ext")  # The thing we're trying to solve for
# d_start = sympy.symbols("d_start")
# s_last = sympy.symbols("s_last")

# vs = v0 * (1.0 - s) + v1 * s  # The velocity that needs cancelling
 
def P(dist):
   return (inner - outer) / (dist - inner) - dist / (inner - outer)
 
def F(dist):
   return (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner))
 
impulse_to_zero = Meff * (vs + F_ext/Meff*delta_s*dt)
impulse_to_zero = Meff * (v0 + s*(vs - v0)/s_last)

# ds = d0 - vs * s * dt
delta_impulse = P(ds+vs*delta_s*dt) - P(d0)
# impulse = P(ds)

impulse = P(d0) - P(d0 + s*(ds - d0)/s_last)

print("Impulse")
print(impulse)

print("Delta Impulse")
print(delta_impulse)

print("Delta Impulse grad")
print(sympy.diff(delta_impulse, delta_s))

eq = impulse - impulse_to_zero
 
print("eq")
print(eq)
 
diff = sympy.diff(eq, s)
 
print("diff")
print(diff)
 
print("Newton")
print(-eq / diff)
 
print("ds")
print(ds)
 
v1 = -F(ds) * (1 - s) * dt / Meff  # TODO Can this be calculated implicitly so we're not using the initial force all the way through the timestep?
 
d1 = ds + 0.5 * -F(ds) / Meff * ((1 - s) * dt)**2  # TODO same as for calculating v1
 
d_d = d1 - d0
d_v = v1 - v0
 
F_update = d_v * Meff / dt
 
print("F update")
print(F_update)
 
print("Impulse update")
print(F_update / dt)
 
d_linear = d0 - d_v * dt
 
F_correction = (d_linear - d1) / dt / dt / Meff
print("F correction")
print(F_correction)
 
"""print(F(d0 - vs * s * dt))
# sF = sympy.integrate(F(d0 - vs * s * dt), s)
sF = 1.0*s/(1.0*inner - 1.0*outer) + 0.25*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1)*log(1.0*s + (-64.0*d0**2*v0**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 128.0*d0**2*v0*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 64.0*d0**2*v1**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 32.0*d0*dt*v0**3*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 32.0*d0*dt*v0**2*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 128.0*d0*inner*v0**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 256.0*d0*inner*v0*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 128.0*d0*inner*v1**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 4.0*dt**2*v0**4*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 32.0*dt*inner*v0**3*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 32.0*dt*inner*v0**2*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 64.0*inner**2*v0**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 128.0*inner**2*v0*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 64.0*inner**2*v1**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 32.0*inner*v0**2 + 32.0*inner*v0*v1 + 32.0*outer*v0**2 - 32.0*outer*v0*v1)/(64.0*inner*v0**2 - 128.0*inner*v0*v1 + 64.0*inner*v1**2 - 64.0*outer*v0**2 + 128.0*outer*v0*v1 - 64.0*outer*v1**2)) - 0.25*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1)*log(1.0*s + (64.0*d0**2*v0**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 128.0*d0**2*v0*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 64.0*d0**2*v1**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 32.0*d0*dt*v0**3*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 32.0*d0*dt*v0**2*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 128.0*d0*inner*v0**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 256.0*d0*inner*v0*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 128.0*d0*inner*v1**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 4.0*dt**2*v0**4*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 32.0*dt*inner*v0**3*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 32.0*dt*inner*v0**2*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 64.0*inner**2*v0**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 128.0*inner**2*v0*v1*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) + 64.0*inner**2*v1**2*sqrt(1/(dt*(-1.0*d0*v0 + 1.0*d0*v1 + 0.25*dt*v0**2 + 1.0*inner*v0 - 1.0*inner*v1)**3))*(1.0*inner - 1.0*outer)*(1.0*v0 - 1.0*v1) - 32.0*inner*v0**2 + 32.0*inner*v0*v1 + 32.0*outer*v0**2 - 32.0*outer*v0*v1)/(64.0*inner*v0**2 - 128.0*inner*v0*v1 + 64.0*inner*v1**2 - 64.0*outer*v0**2 + 128.0*outer*v0*v1 - 64.0*outer*v1**2)) + (-1.0*inner*v0 + 1.0*outer*v0 + s*(2.0*inner*v0 - 2.0*inner*v1 - 2.0*outer*v0 + 2.0*outer*v1))/(-4.0*d0**2*v0 + 4.0*d0**2*v1 + 1.0*d0*dt*v0**2 + 8.0*d0*inner*v0 - 8.0*d0*inner*v1 - 1.0*dt*inner*v0**2 - 4.0*inner**2*v0 + 4.0*inner**2*v1 + s**2*(-4.0*d0*dt*v0**2 + 8.0*d0*dt*v0*v1 - 4.0*d0*dt*v1**2 + 1.0*dt**2*v0**3 - 1.0*dt**2*v0**2*v1 + 4.0*dt*inner*v0**2 - 8.0*dt*inner*v0*v1 + 4.0*dt*inner*v1**2) + s*(4.0*d0*dt*v0**2 - 4.0*d0*dt*v0*v1 - 1.0*dt**2*v0**3 - 4.0*dt*inner*v0**2 + 4.0*dt*inner*v0*v1))
print("sF:")
print(sF)
 
vel_delta = (sF - sF.subs(s, 0)) / Meff
 
eq = vel_delta - vs
 
print("eq:")
print(eq)
 
print("simp")
print(sympy.simplify(eq))
 
diff = (F(d0 - vs * s * dt) - F(d0 - vs * s * dt).subs(s, 0)) / Meff - sympy.diff(vs, s)
 
print("diff")
print(diff)"""

