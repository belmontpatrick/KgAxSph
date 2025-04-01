import numpy as np
from parameters import LX,LR
import scipy.special as sp

# Basis function in r using Odd Chebyshev Polynomials

def SB(n, r):
    res = np.sin((2*n+1)*np.arctan(LR/r)) 
    return res

def dSB(n,r):
    res = -np.cos((2*n+1)*np.arctan(LR/r))*(2*n+1)*LR/(r**2*(1+LR**2/r**2)) 
    return res


def ddSB(n,r):
    res = (-np.sin((2*n+1)*np.arctan(LR/r))*(2*n+1)**2*LR**2/(r**4*(1+LR**2/r**2)**2)+
2*np.cos((2*n+1)*np.arctan(LR/r))*(2*n+1)*LR/(r**3*(1+LR**2/r**2))-2*np.cos((2*n+1)*np.arctan(LR/r))*(2*n+1)*LR**3/(r**5*(1+LR**2/r**2)**2))
    return res

# Basis function in x using Legendre Polynomials

def P(i, x):
    return sp.legendre(i)(x)
# print(P(12,5))

def dP(i, x):
    return sp.legendre(i).deriv()(x)

def ddP(i, x):
    return sp.legendre(i).deriv().deriv()(x)
