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

# rcol = np.loadtxt("rcolPR60L5.txt", delimiter="\t")

# " Base Matrix (Tchebyshev Polinomials) in r: "

# SB_r = np.zeros([PR+1,PR+1])
# dSB_r = np.zeros([PR+1,PR+1])
# ddSB_r = np.zeros([PR+1,PR+1])

# for i in range(PR+1):
#   SB_r[i,] = SB(i,rcol)                                                 

# for i in range(PR+1):
#   dSB_r[i,] = dSB(i,rcol)

# for i in range(PR+1):
#   ddSB_r[i,] = ddSB(i,rcol)

"Importing Collocation Points in x:"

xcol = np.loadtxt("xcolPX10_Henrique.txt", delimiter="\t")

"Base Matrix (Legendre Polinomials) in x: "

P_x = np.zeros([PX+1,PX+1])
dP_x = np.zeros([PX+1,PX+1])
ddP_x = np.zeros([PX+1,PX+1])

for i in range(PX+1):
  P_x[i,] = P(i,xcol)

for i in range(PX+1):
  dP_x[i,] = dP(i,xcol)

for i in range(PX+1):
  ddP_x[i,] = ddP(i,xcol)

print(P_x)