import numpy as np
from scipy.optimize import fsolve
import scipy.special as sp

# Your existing code
PX = 10
P_prime = sp.legendre(2 * PX + 3).deriv()
x_roots = fsolve(P_prime, np.cos(np.pi * (np.arange(1, 2 * PX + 3) / (2 * PX + 3))))
x_col_prel = np.sort(x_roots)
x_col = np.flip(x_roots[:PX + 1])

np.savetxt("xcolPX10.txt", x_col, fmt="%.16e", delimiter="\t")