import numpy as np
# import collocation as cp
import collocation as cm
import parameters as par

r_col = cm.r_col
x_col = cm.x_col
Psi = cm.Psi
Psi_inv = cm.inv_psi
drPsi = cm.rPsi
ddrPsi = cm.rrPsi
dxPsi = cm.xPsi
ddxPsi = cm.xxPsi
N = par.N
px = par.px


r = np.tile(r_col.repeat(px + 1).reshape(-1, 1), (1, (N + 1) * (px + 1)))
x = np.tile(np.tile(x_col, N + 1).reshape(-1, 1), (1, (N + 1) * (px + 1)))

def dda(c):

    # # drphi = np.dot(drPsi,c)
    # # drrphi = np.dot(ddrPsi,c)
    # # dxxphi = np.dot(ddxPsi,c)
    # # dxphi = np.dot(dxPsi,c)

    # RHS = ddrPsi + 2/r*drPsi  + (-x**2 + 1) * ddxPsi/r**2 - 2*x * dxPsi

    # # np.dot(np.dot(c, col.ddrPsi + 2/r*col.drPsi + (-x**2 + 1)*col.ddxPsi/r**2 - 2*x*col.dxPsi), col.Psi)

    # return np.dot(RHS,Psi_inv)
    return np.dot(np.dot(c, ddrPsi + 2/r*drPsi + (-x**2 + 1)*ddxPsi/r**2 - 2*x*dxPsi), Psi_inv)

def phi(c):
    return np.dot(Psi,c)