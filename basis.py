import time
import numpy as np
from parameters import L0
import scipy.special as sp

##BASES
#r basis
def SB(n, r):
    return np.sin((n+1)*np.arctan(L0/r))

# np.sin((n + 1) * ((np.pi / 2) - np.arctan(r / L0)))

# np.sin((n+1)*np.arctan(L0/r))
# print(SB(1,2))

def rSB(n, r):
    return -np.cos((n+1)*np.arctan(L0/r))*(n+1)*L0/(r**2*(1+L0**2/r**2)) 

# - np.cos((n + 1) * ((np.pi / 2) - np.arctan(r / L0))
#                    ) * (n + 1) / (L0 * (1 + (r / L0)**2))

# -np.cos((n+1)*np.arctan(L0/r))*(n+1)*L0/(r**2*(1+L0**2/r**2)) 

def rrSB(n, r):
    return (-np.sin((n+1)*np.arctan(L0/r))*(n+1)**2*L0**2/(r**4*(1+L0**2/r**2)**2)+2*np.cos((n+1)*np.arctan(L0/r))*(n+1)*L0/(r**3*(1+L0**2/r**2))-2*np.cos((n+1)*np.arctan(L0/r))*(n+1)*L0**3/(r**5*(1+L0**2/r**2)**2))

# (2 * np.cos((n + 1) * (np.pi / 2 - np.arctan(r / L0))) *
#            (n + 1) * r / (L0**3 * (1 + (r / L0)**2)**2)
#            - np.sin((n + 1) * (np.pi / 2 - np.arctan(r / L0))) *
#            (n + 1)**2 / (L0**2 * (1 + (r / L0)**2)**2)
#            )

# (-np.sin((n+1)*np.arctan(L0/r))*(n+1)**2*L0**2/(r**4*(1+L0**2/r**2)**2)+2*np.cos((n+1)*np.arctan(L0/r))*(n+1)*L0/(r**3*(1+L0**2/r**2))-2*np.cos((n+1)*np.arctan(L0/r))*(n+1)*L0**3/(r**5*(1+L0**2/r**2)**2))


# #x basis
def P(i, x):
    return sp.legendre(i)(x)

# Legendre.basis(i)(x) 
# # sp.eval_legendre(i, x)
# # sp.legendre(i)(x)
# # print(P(12,5))

def xP(i, x):
    return sp.legendre(i).deriv()(x)
# poly = Legendre.basis(i).deriv(1) returuno poly(x)
# # sp.eval_legendre(i, x, derivative=1)
# # sp.legendre(i).deriv()(x)

def xxP(i, x):
    return sp.legendre(i).deriv().deriv()(x)