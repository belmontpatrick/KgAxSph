import numpy as np
from parameters import L0, N, px
from basis import SB, rSB, rrSB, P, xP, xxP
import scipy.special as sp
from scipy.optimize import fsolve
from numpy.polynomial.legendre import Legendre

##COLLOCATION POINTS

#new r collocation points with linspace
k_values = np.linspace(0, 2*N + 3, 2*N + 4, dtype=np.float64)
x__col = np.cos(np.pi * k_values / (2*N + 3))
epsilon = 1e-15  # Valor pequeno para evitar divisão por zero
r_col_pre = L0 * x__col / np.sqrt(1 - x__col**2 + epsilon)
r_col = np.flip(np.array([r_col_pre[N + 2 - k] for k in range(1, N + 2)])) #ordem inversa com np.flip

#collocation points for x
# P_prime = sp.legendre(2 * px + 3).deriv()
# x_roots = fsolve(P_prime, np.cos(np.pi * (np.arange(1, 2 * px + 3) / (2 * px + 3))))
# x_col_prel = np.sort(x_roots)
# x_col = -x_col_prel[:px + 1]

gegen_roots, _ = sp.roots_gegenbauer( 2 * px + 3 - 1,3/2)
gegen_col_prel = np.sort(gegen_roots)
x_col = - gegen_col_prel[:px + 1]

#collocation points on the bases
# r basis
SB_ = np.zeros([N+1,N+1])
SB_r = np.zeros([N+1,N+1])
SB_rr = np.zeros([N+1,N+1])

SB_0 = np.zeros([1,N+1])

for i in range(N+1):
    SB_[i,] = SB(2 * i,r_col)
# np.savetxt('SB_in', SB_.T, fmt='%.20f') #confere

for i in range(N+1):
    SB_r[i,] = rSB(2 * i,r_col)
# np.savetxt('SBr_in', SB_r, fmt='%.20f') #confere

for i in range(N+1):
    SB_rr[i,] = rrSB(2 * i,r_col)
# np.savetxt('SBrr_in', SB_rr, fmt='%.20f') #confere

# for i in range(N+1):
#     SB_0[i,] = SB(i,1)
    
# print(SB_0)

#x basis
P_ = np.zeros((px + 1, px + 1))
P_x = np.zeros((px + 1, px + 1))
P_xx = np.zeros((px + 1, px + 1))

for i in range(px + 1):
    P_[i,] = P(2*i,x_col)
# np.savetxt('P_in', P_, fmt='%.20f') #confere
    
for i in range(px + 1):
    P_x[i,] = xP(2*i,x_col)
# np.savetxt('Px_in', P_x, fmt='%.20f') #confere
    
for i in range(px + 1):
    P_xx[i,] = xxP(2*i,x_col)

Psi    = np.kron(SB_.T,P_.T)
rPsi   = np.kron(SB_r.T,P_.T)  #/L
rrPsi  = np.kron(SB_rr.T,P_.T)  #/L**2
xPsi   = np.kron(SB_.T,P_x.T)
xxPsi  = np.kron(SB_.T,P_xx.T)

inv_psi = np.linalg.inv(Psi)