import numpy as np
import basis as bs
import parameters as par
import matplotlib.pyplot as plt

PR = par.PR
PX = par.PX

a_datatest = np.loadtxt("resultados_rk4_tst.txt")

atest = a_datatest[99, 1:]

rplot = np.linspace(0.000001, 10, 100)
xplot = np.linspace(0.000001, 10, 100)

Rplot, Xplot = np.meshgrid(rplot, xplot, indexing='ij')

def phi_approx(c0,r,x):
    res = sum(c0[k] * bs.Basis(k, r, x) for k in range((PR+1) * (PX+1)))
    return res

phi50 = phi_approx(Rplot, Xplot, a_datatest)

# Assuming phi50, Rplot, and Zplot are already defined
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create the 3D surface plot
surf = ax.plot_surface(Rplot, Xplot, phi50, cmap='viridis')

# Add labels and a color bar
ax.set_xlabel('R')
ax.set_ylabel('Z')
ax.set_zlabel('Phi')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)

# Show the plot
plt.show()