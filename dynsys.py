import numpy as np
# import collocation as cp
import collocation_matrices as cm
import parameters as par

rcol = cm.rcol
xcol = cm.xcol
Psi = cm.Psi
Psi_inv = cm.Psi_inv
drPsi = cm.drPsi
ddrPsi = cm.ddrPsi
dxPsi = cm.dxPsi
ddxPsi = cm.ddxPsi
PR = par.PR
PX = par.PX

# Rm = np.repeat(rcol,len(xcol))
# Rr = Rm.reshape(-1,1)

r = np.tile(rcol.repeat(PX + 1).reshape(-1, 1), (1, (PR + 1) * (PX + 1)))

x = np.tile(xcol.repeat(PX + 1).reshape(-1, 1), (1, (PR + 1) * (PX + 1)))

# Xm = np.repeat(xcol,len(rcol))
# Xr = Xm.reshape(-1,1)

def dda(c):

    drphi = np.dot(drPsi,c)
    drrphi = np.dot(ddrPsi,c)
    dxxphi = np.dot(ddxPsi,c)
    dxphi = np.dot(dxPsi,c)

    RHS = drrphi + drphi * 2 /r  + (1/ ( r ) ** 2) * ((1 - x ** 2) * dxxphi * 2 * x * dxphi)

    return np.dot(Psi_inv,RHS)

def phi(c):
    return np.dot(Psi,c)
