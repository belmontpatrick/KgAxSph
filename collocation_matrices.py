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

"Importing Collocation Points in r:"

rcol = np.loadtxt("rcolPR60L1.txt", delimiter="\t")

" Base Matrix (Tchebyshev Polinomials) in r: "

SB_r = np.zeros([PR+1,PR+1])
dSB_r = np.zeros([PR+1,PR+1])
ddSB_r = np.zeros([PR+1,PR+1])

for i in range(PR+1):
  SB_r[i,] = SB(i,rcol)                                                 

for i in range(PR+1):
  dSB_r[i,] = dSB(i,rcol)

for i in range(PR+1):
  ddSB_r[i,] = ddSB(i,rcol)

"Importing Collocation Points in x:"

xcol = np.loadtxt("xcolPX30.txt", delimiter="\t")

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

### Multiplication of the Matrices:

# Psi = np.tile(SB_r.T,(PX+1,PX+1)) * repelem(P_x.T,(PR+1,PR+1))

Psi = np.kron(SB_r.T, P_x.T)

Psi_inv = np.linalg.inv(Psi)

# np.savetxt("PsiPR60L1PX30.txt", Psi, fmt="%.16e", delimiter="\t")

drPsi = np.kron(dSB_r.T, P_x.T)

ddrPsi = np.kron(ddSB_r.T, P_x.T)

dxPsi = np.kron(SB_r.T, dP_x.T)

ddxPsi = np.kron(SB_r.T, ddP_x.T)