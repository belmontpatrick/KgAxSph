import numpy as np
from initialdata import a0, da
# from collocation import Psi, inv_psi, rPsi, rrPsi, xPsi, xxPsi,r_col, x_col
from timegrid import h, It, t
from parameters import N, px
from dynsys import dda, phi

# r = np.tile(r_col.repeat(px + 1).reshape(-1, 1), (1, (N + 1) * (px + 1)))
# x = np.tile(np.tile(x_col, N + 1).reshape(-1, 1), (1, (N + 1) * (N + 1)))

with open('resultados_a_def_tst.txt', 'w') as a0_file, \
     open('resultados_da_def_tst.txt', 'w') as da_file:

    for i in range(It):  # Runge Kutta 4th order

        phi_ = phi(a0)
        dda_ = dda(a0)
        L1 = h*(da)
        K1 = h*(dda_)

        phi_ = phi(a0 + L1/2)
        dda_ = dda(a0 + L1/2)
        L2 = h*(da + K1/2)
        K2 = h*(dda_)

        phi_ = phi(a0 + L2/2)
        dda_ = dda(a0 + L2/2)
        L3 = h*(da + K2/2)
        K3 = h*(dda_)

        phi_ = phi(a0 + L3)
        dda_ = dda(a0 + L3)
        L4 = h*(da + K3)
        K4 = h*(dda_)

        da = da + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
        a0 = a0 + 1/6 * (L1 + 2*L2 + 2*L3 + L4)
        
        a0_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in a0.flatten()])
        da_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in da.flatten()])
        
        a0_file.write(a0_line + "\n")
        da_file.write(da_line + "\n")

# r = np.tile(rcol.repeat(PX + 1).reshape(-1, 1), (1, (PR + 1) * (PX + 1)))
# x = np.tile(np.tile(xcol, PR + 1).reshape(-1, 1), (1, (PR + 1) * (PX + 1)))


# with open('resultados_a0_gegen.txt', 'w') as a0_file, \
#      open('resultados_da_gegen.txt', 'w') as da_file:

#     for i in range(It):  # Runge Kutta 4th order

#         phi = np.dot(a0, Psi)
#         dda = np.dot(np.dot(a0, rrPsi + 2/r*rPsi + (-x**2 + 1)*xxPsi/r**2 - 2*x*xPsi), inv_psi)
#         L1 = h*(da)
#         K1 = h*(dda)

#         phi = np.dot(a0 + L1/2, Psi)
#         dda = np.dot(np.dot(a0 + L1/2, rrPsi + 2/r*rPsi + (-x**2 + 1)*xxPsi/r**2 - 2*x*xPsi), inv_psi)
#         L2 = h*(da + K1/2)
#         K2 = h*(dda)

#         phi = np.dot(a0 + L2/2, Psi)
#         dda = np.dot(np.dot(a0 + L2/2, rrPsi + 2/r*rPsi + (-x**2 + 1)*xxPsi/r**2 - 2*x*xPsi), inv_psi)
#         L3 = h*(da + K2/2)
#         K3 = h*(dda)

#         phi = np.dot(a0 + L3, Psi)
#         dda = np.dot(np.dot(a0 + L3, rrPsi + 2/r*rPsi + (-x**2 + 1)*xxPsi/r**2 - 2*x*xPsi), inv_psi)
#         L4 = h*(da + K3)
#         K4 = h*(dda)

#         da = da + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
#         a0 = a0 + 1/6 * (L1 + 2*L2 + 2*L3 + L4)
        
#         a0_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in a0.flatten()])
#         da_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in da.flatten()])
        
#         a0_file.write(a0_line + "\n")
#         da_file.write(da_line + "\n")

