import numpy as np
import scipy.special as sp
import parameters as par
import basis as bs

PR = par.PR
LR = par.LR
PX = par.PX
P = bs.P
dP = bs.dP
ddP = bs.ddP
SB = bs.SB
dSB = bs.dSB
ddSB = bs.ddSB

"Legendre Collocation Points in x:"

def P_colpoints(n):
  "getting gegenbauer roots"
  gegen_roots, _ = sp.roots_gegenbauer( 2 * n + 3 - 1,3/2)
  gegen_col_prel = np.sort(gegen_roots)
  return - gegen_col_prel[:n + 1]

xcol = P_colpoints(PX)

"Base Matrix (Legendre Polinomials) in x: "

P_x = np.zeros([PX+1,PX+1])
dP_x = np.zeros([PX+1,PX+1])
ddP_x = np.zeros([PX+1,PX+1])

for i in range(PX+1):
  P_x[i,] = P(2*i,xcol)

for i in range(PX+1):
  dP_x[i,] = dP(2*i,xcol)

for i in range(PX+1):
  ddP_x[i,] = ddP(2*i,xcol)

P_x_inv = np.linalg.inv(P_x)

"Chebyshev Collocation Points in r:"

def SB_colpoints(n):
  k_values = np.linspace(0, 2*PR + 3, 2*PR + 4, dtype=np.float64)
  x__col = np.cos(np.pi * k_values / (2*PR + 3))
  epsilon = 1e-15  # Valor pequeno para evitar divisão por zero
  r_col_pre = LR * x__col / np.sqrt(1 - x__col**2 + epsilon)
  return np.flip(np.array([r_col_pre[PR + 2 - k] for k in range(1, PR + 2)]))

rcol = SB_colpoints(PR)

" Base Matrix (Tchebyshev Polinomials) in r: "

SB_r = np.zeros([PR+1,PR+1])
dSB_r = np.zeros([PR+1,PR+1])
ddSB_r = np.zeros([PR+1,PR+1])

for i in range(PR+1):
  SB_r[i,] = SB(2 * i,rcol)                                                 

for i in range(PR+1):
  dSB_r[i,] = dSB(2 * i,rcol)

for i in range(PR+1):
  ddSB_r[i,] = ddSB(2 * i,rcol)

" Big Matrix: "

Psi = np.kron(SB_r.T, P_x.T)

Psi_inv = np.linalg.inv(Psi)

drPsi = np.kron(dSB_r.T, P_x)

ddrPsi = np.kron(ddSB_r.T, P_x.T)

dxPsi = np.kron(SB_r.T, dP_x.T)

ddxPsi = np.kron(SB_r.T, ddP_x.T)

# np.savetxt("PsiPR20L1PX20.txt", Psi, fmt="%.16e", delimiter="\t")
# print("Psi saved ")

# np.savetxt("Psi_invPR20L1PX20.txt", Psi_inv, fmt="%.16e", delimiter="\t")
# print("Psi_inv saved ")