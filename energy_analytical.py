import sympy as sp
import numpy as np
from parameters import A0

# Define symbols
r, t = sp.symbols('r t', real=True, positive=True)
x = sp.symbols('x', real=True)

# Field ansatz: phi = A0 * exp(-r^2) * (1 - cos(x))
phi = A0 * sp.exp(-r**2) * (1 - sp.cos(x))

# Time derivative
dt_phi = 0

# Spatial derivatives
dr_phi = sp.diff(phi, r)      # ∂ϕ/∂r
dx_phi = sp.diff(phi, x)      # ∂ϕ/∂x

# Hamiltonian density integrand
integrand = r**2 * (sp.Abs(dt_phi)**2 + r**2 * (dr_phi)**2 + (1 - x**2) * (dx_phi))**2

# First integrate over x from -1 to 1
integrand_r = sp.integrate(integrand, (x, -1, 1))

# Now integrate over r from 0 to infinity
H = sp.Integral(sp.pi * integrand_r, (r, 0, sp.oo))

# Compute the analytical result
analytical_energy = H.doit()

# Numerical evaluation (optional)
eval_energy = sp.N(analytical_energy)

#print("Analytical Hamiltonian:")
#sp.pprint(analytical_energy)
print("Evaluated  analytical energy:")
print(eval_energy)