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
    # return A * np.exp(-r**2)*(1 - x ** 2)
    # return A * np.exp(-r**2 - x**2)
   return A * np.exp(-r**2) * (1 - (np.cos(x))**2)

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

r_plot = np.linspace(0.0001, 5, 100) 
x_plot = np.linspace(-0.5, 0.5, 100)

Rplot, Xplot = np.meshgrid(r_plot,x_plot)

da0sph = np.zeros((PR + 1) * (PX + 1))

error = gaussian_exact(A0,Rplot,Xplot) - phi_approx(a0sph,Rplot,Xplot)

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

fig3 = plt.figure()
ax3 = plt.axes(projection = '3d')
ax3.plot_surface(Rplot, Xplot, error, cmap='viridis')
ax3.set_xlabel('r')
ax3.set_ylabel('x')
ax3.set_zlabel('error')
plt.show()

####
