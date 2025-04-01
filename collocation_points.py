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

# def P_colpoints(n):
#     P_prime = sp.legendre(2 * n + 1).deriv()
#     x_roots = fsolve(P_prime, np.cos(np.pi * (np.arange(1, 2 * n + 2) / (2 * n + 2))))
#     x_col_prel = np.concatenate([[-1], np.sort(x_roots), [1]])
#     return -np.flip(x_col_prel[:PX + 1])

# P_prime = sp.legendre(2 * PX + 1).deriv()
# x_roots = fsolve(P_prime, np.cos(np.pi * (np.arange(1, 2 * PX + 2) / (2 * PX + 2))))
# x_col_prel = np.concatenate(([-1], np.sort(x_roots), [1]))
# xcol = -np.flip(x_col_prel[:PX + 1])

xcol_Henrique = np.array([.668379937372285781e-1,
.199321253390832667,
.328247613375510912,
.451316373214322618,
.566331357979295312,
.671240105264128700,
.764170482420493308,
.843464070154872041,
.907705675113506522,
.955748220929886358,
.986730553505160884,
])

np.savetxt("xcolPX10_Henrique.txt", xcol_Henrique, fmt="%.16e", delimiter="\t")

print("xcol_Henrique", xcol_Henrique)
def SB_colpoints(n):
    return np.cos(np.arange(2*n + 4)*np.pi /(2*n + 3))

"Chebyshev collocation points in r"

col_r = SB_colpoints(PR)

colr = col_r[1:PR+2]

colr1 = LR * colr/(np.sqrt(1-colr**2))                       # physical domain 
rcol = np.flip(colr1)

# np.savetxt("rcolPR60L5.txt", rcol, fmt="%.16e", delimiter="\t")