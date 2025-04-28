import numpy as np
import parameters as par
import basis as bs
import matplotlib.pyplot as plt
from matplotlib import cm
import collocation_matrices as col

A0 = par.A0
rcol = col.rcol
xcol = col.xcol
PX = par.PX
PR = par.PR
Basis = bs.Basis
P = bs.P
SB = bs.SB
Psi = col.Psi
Psi_inv = col.Psi_inv

"Defining initial data"

def gaussian_exact(A,r,x):
    return A * np.exp(-r**2)*(1 - x ** 2)

phi_exact_col_arr = [gaussian_exact(A0,col.rcol[j], col.xcol[k])
                          for j in range(PR + 1)
                          for k in range(PX + 1)]

#a0_col_arr = np.dot(Psi_inv, phi_exact_col_arr)

a0sph = np.dot(Psi_inv, phi_exact_col_arr)

#np.savetxt("a0sphPR3L1PX3.txt", a0sph, fmt="%.16e", delimiter="\t")

def phi_approx(c0,r,x):
    res = sum(c0[k] * Basis(k, r, x) for k in range((PR+1) * (PX+1)))
    return res

#### PLOTS #####

r_plot = np.linspace(0, 10, 100) 
x_plot = np.linspace(-1, 1, 100)

Rplot, Xplot = np.meshgrid(r_plot,x_plot)

da0sph = np.zeros(((PR + 1) * (PX + 1),1))


Z = gaussian_exact(A0,Rplot,Xplot)
# Y = phi_approx(a0_col_arr,Rplot,Xplot)
Y = phi_approx(a0sph,Rplot,Xplot)

fig1 = plt.figure()
ax1 = plt.axes(projection = '3d')
ax1.plot_surface(Rplot, Xplot, Z, cmap='viridis')
ax1.set_xlabel('r')
ax1.set_ylabel('x')
ax1.set_zlabel('phi0_exact')
# plt.savefig('exact_phi_trunc121.png', dpi=300, bbox_inches='tight')

fig2 = plt.figure()
ax2 = plt.axes(projection = '3d')
ax2.plot_surface(Rplot, Xplot, Y, cmap='viridis')
ax2.set_xlabel('r')
ax2.set_ylabel('x')
ax2.set_zlabel('phi0_approx')
# plt.savefig('approx_phi_trunc121.png', dpi=300, bbox_inches='tight')

# plt.figure(figsize=(5, 5))
# fig3 = plt.figure()
# ax3 = plt.axes(projection = '3d')
# ax3.plot_surface()
plt.show()

####

# Pcounts = (PR + 1) * (PX + 1)

# ones = np.ones((PX+1))

# P_1 = np.zeros([PX+1,PX+1])

# for i in range(PX+1):
#   P_1[i,] = P(2*i,ones)

# zeros = np.zeros(rcol.shape)


# SB_0 = np.zeros([PR+1,PR+1])

# for i in range(PR+1):
#   SB_0[i,] = SB(i,zeros)

# print(SB_0.T[0])

# print(P_1.T[0])

# basis01 = np.kron(SB_0.T, P_1.T)

# b01 = basis01[0]

# # print(b01)

# phi001 = sum(a0sph[k] * b01[k] for k in range((PR+1) * (PX+1)))

# print(phi001)

# rplot = np.linspace(0.00001, 1, 20)
# xplot = np.linspace(-1, 1, 20)

# Rplot, Xplot = np.meshgrid(rplot, xplot)

# #phi0_approx = phi_approx(Rplot, Xplot, a0sph)
# phi0_approx = phi_approx(a0_col_arr, Rplot, Xplot)

# phi0_approx =  [alpha_m1_exact(col.col_r_shift[j], col.col_z_shift[k], A0_TEST, SIGMA_R, SIGMA_Z)
#                          for j in range(PR_TEST + 1) for k in range(PZ_TEST + 1)]

# phiinidata = iniaxsph(A0, Rplot, np.flip(Xplot))

# ax = plt.axes(projection = '3d')
# # ax.plot_surface(Rplot, Xplot, phiinidata, cmap=cm.viridis)
# ax.plot_surface(Rplot, Xplot, phiinidata, cmap=cm.viridis)
# plt.show()

# fig1 = plt.figure()
# ax1 = plt.axes(projection="3d")
# R, Z = np.meshgrid(np.linspace(0.0001, 10, 10), np.linspace(0.0001, 10, 10))
# Y1 = iniaxsph(A0,Rplot, Xplot)
# ax1.plot_surface(R, Z, Y1, cmap='plasma')
# ax1.set_xlabel('r')
# ax1.set_ylabel('z')
# ax1.set_zlabel('alpha0_exact')

# fig2 = plt.figure()
# ax2 = plt.axes(projection="3d")
# R, Z = np.meshgrid(np.linspace(0.0001, 10,0), np.linspace(0.0001, 10, 10))
# Y2 = phi_approx(R, Z, phi0)
# ax2.plot_surface(R, Z, Y2, cmap='plasma')
# ax2.set_xlabel('r')
# ax2.set_ylabel('z')
# ax2.set_zlabel('alpha0_approx')