import numpy as np
import basis as bs
import parameters as par
import scipy.special as sp
from scipy.optimize import fsolve

PR = par.PR
LR = par.LR
PX = par.PX
P = bs.P
dP = bs.dP
ddP = bs.ddP
SB = bs.SB
dSB = bs.dSB
ddSB = bs.ddSB


"Legendre Collocation points in x"

def P_colpoints(n):
    P_prime = sp.legendre(2 * PX + 3).deriv()
    x_roots = fsolve(P_prime, np.cos(np.pi * (np.arange(1, 2 * PX + 3) / (2 * PX + 3))))
    x_col_prel = np.sort(x_roots)
    x_col = np.flip(x_roots[:PX + 1])

    return x_col

xcol = P_colpoints(PX)

np.savetxt("xcolPX30.txt", xcol, fmt="%.16e", delimiter="\t")

# print("xcol_Henrique", xcol_Henrique)
def SB_colpoints(n):
    return np.cos(np.arange(2*n + 4)*np.pi /(2*n + 3))

"Chebyshev collocation points in r"

col_r = SB_colpoints(PR)

colr = col_r[1:PR+2]

colr1 = LR * colr/(np.sqrt(1-colr**2))                       # physical domain 
rcol = np.flip(colr1)

np.savetxt("rcolPR60L1.txt", rcol, fmt="%.16e", delimiter="\t")