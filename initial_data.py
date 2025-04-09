import collocation_matrices as cm
import numpy as np
import parameters as par

A0 = par.A0
Psi_inv = cm.Psi_inv
P_x = cm.P_x
PX = par.PX
PR = par.PR

"Importing Collocation Points:"
rcol = np.loadtxt("rcolPR60L1.txt", delimiter="\t")

xcol = np.loadtxt("xcolPX30.txt", delimiter="\t")

rtile = np.tile(rcol, (len(xcol)))

xrepeat = np.repeat(xcol, len(rcol))

"Defining inital data"

def iniaxsph(A,r,x):
    return A * np.exp(-r**2)*(1 - (np.cos(x)) ** 2)

inidata = iniaxsph(A0, rtile, xrepeat)

a0sph = np.dot(Psi_inv, inidata)

P = (PR+1) * (PX + 1)

da0sph = np.zeros(P)

print(a0sph)