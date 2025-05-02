import numpy as np
import matplotlib.pyplot as plt
# import math
import scipy.special as sp
import scipy
from scipy.optimize import fsolve
from numpy.polynomial.legendre import Legendre

np.set_printoptions(precision=16)

#Parameters
N = 4
L0 = 1
# SIGMA_r = 1
A0 = 0.002
# r0 = 2
px = 2