import numpy as np
import basis as bs
import parameters as par
import matplotlib.pyplot as plt
import timegrid as gr

PR = par.PR
PX = par.PX

a_datatest = np.loadtxt("resultados_a_direct20_20.txt")

a_test = a_datatest[0, 1:]
# print(a_TRUNC_2_2.shape)

print(a_test.shape)

rplot = np.linspace(0.000001, 5, 100)
xplot = np.linspace(-1, 1, 100)

Rplot, Xplot = np.meshgrid(rplot, xplot, indexing='ij')

def phi_approx(c0,r,x):
    res = sum(c0[k] * bs.Basis(k, r, x) for k in range((PR+1) * (PX+1)))
    return res


# phi_test = np.zeros(gr.It, len(rplot), len(xplot))

# for i in range(len(a_test)):
#     phi_test = phi_approx(a_test[i], Rplot, Xplot)

phitest = phi_approx(a_test,Rplot, Xplot)

# Assuming phi50, Rplot, and Zplot are already defined
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create the 3D surface plot
surf = ax.plot_surface(Rplot, Xplot, phitest, cmap='viridis')

# Add labels and a color bar
ax.set_xlabel('R')
ax.set_ylabel('X')
ax.set_zlabel('Phi')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)

# Show the plot
plt.show()