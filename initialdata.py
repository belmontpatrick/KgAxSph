import numpy as np
from collocation import Psi, inv_psi, r_col,x_col
from parameters import L0, N, px, A0
import time
from basis import SB, P

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Rt, Xt = np.meshgrid(r_col, x_col)
# R = np.tile(r_col,(len(x_col)))
# X = np.repeat(x_col,(len(r_col)))

#Definir os intervalos
r_plot = np.linspace(0, 5, 100) 
x_plot = np.linspace(-0.5, 0.5, 100)

Rplot, Xplot = np.meshgrid(r_plot,x_plot)

# def gaussian(r,x):
#     return A0 * np.exp(-r**2) * (-x**2 + 1)

def gaussian_exact(A0,r,x):
    return A0 * np.exp(-r**2) * (1 - (np.cos(x))**2)

# A0 * np.exp(-r**2 -x**2)
# A0 * np.exp(-r**2) * (1 - x**2)


# gaussian_col = gaussian_exact(A0,R,X)
# gaussian_new = gaussian_col.reshape(-1,1)


coef_col = [gaussian_exact(A0, r_col[j], x_col[k])
                         for j in range(N + 1) for k in range(px + 1)]

a0 = np.dot(inv_psi,coef_col)

da = np.zeros((N+1)*(px+1))

def Psi_plot(n, r, x):
    list = [SB(2 * i, r) * P(2 * j, x)
            for i in range(N + 1) for j in range(px + 1)]
    return list[n]

def gaussian_approx(r,x,b0):
    res = sum(b0[k] * Psi_plot(k,r,x) for k in range((N+1) * (px+1)))
    return res

" Plotting the results "

Z = gaussian_exact(A0,Rplot,Xplot)
Y = gaussian_approx(Rplot,Xplot,a0)
# print(Y)
# np.savetxt('gaussian_exact', Z, fmt='%.20f')
# np.savetxt('gaussian_approx', Y, fmt='%.20f')
W = Z - Y
# print(W)


### PLOTS #####

# fig1 = plt.figure()
# ax1 = plt.axes(projection = '3d')
# ax1.plot_surface(Rplot, Xplot, Z, cmap='viridis')
# ax1.set_xlabel('r')
# ax1.set_ylabel('x')
# ax1.set_zlabel('phi0_exact')
# # plt.savefig('exact_phi_trunc121.png', dpi=300, bbox_inches='tight')

# fig2 = plt.figure()
# ax2 = plt.axes(projection = '3d')
# ax2.plot_surface(Rplot, Xplot, Y, cmap='viridis')
# ax2.set_xlabel('r')
# ax2.set_ylabel('x')
# ax2.set_zlabel('phi0_approx')
# # plt.savefig('approx_phi_trunc121.png', dpi=300, bbox_inches='tight')

# fig3 = plt.figure()
# ax3 = plt.axes(projection = '3d')
# ax3.plot_surface(Rplot, Xplot, W, cmap='viridis')
# ax3.set_xlabel('r')
# ax3.set_ylabel('x')
# ax3.set_zlabel('erro')
# # plt.savefig('exact_phi_trunc121.png', dpi=300, bbox_inches='tight')

# plt.show()

