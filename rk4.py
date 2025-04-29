import numpy as np
import basis as bs
import initial_data as ini
import timegrid as gr
# import dynsys as ds
import parameters as par
import collocation_matrices as col

np.set_printoptions(precision=16)

PR = par.PR
PX = par.PX
a0sph = ini.a0sph
da0sph = ini.da0sph
h = gr.h
#t = gr.t
N = gr.N
# phi = ds.phi
# dda = ds.dda
t = gr.t
# r = ds.r
# x = ds.x
rcol = col.rcol
xcol = col.xcol
drPsi = col.drPsi
ddrPsi = col.ddrPsi
dxPsi = col.dxPsi
ddxPsi = col.ddxPsi
Psi = col.Psi
Psi_inv = col.Psi_inv

#Runge-Kutta 4th order


# file_asph = open("rk4_asph.txt", "w")
# file_dacyl = open("rk4_dacyl.txt", "w")
# with open('results_rk4_test_TRUNC2_2.txt', 'w') as file:

# # t = 0.0           # Initial time

#   for i in range(N):  # Runge Kutta 4th order

    # phi_ = phi(a0sph)
    # dda_ = dda(a0sph)
    # L1 = h*da0sph
    # K1 = h*(dda_)

    # phi_ = phi(a0sph + L1/2)
    # dda_ = dda(a0sph + L1/2)
    # L2 = h*(da0sph + K1/2)
    # K2 = h*(dda_)

    # phi_ = phi(a0sph + L2/2)
    # dda_ = dda(a0sph + L2/2)
    # L3 = h*(da0sph + K2/2)
    # K3 = h*(dda_)

    # phi_ = phi(a0sph + L3)
    # dda_ = dda(a0sph + L3)
    # L4 = h*(da0sph + K3)
    # K4 = h*(dda_)

#     da0sph = da0sph + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
#     a0sph = a0sph + 1/6 * (L1 + 2*L2 + 2*L3 + L4)
    
#     line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in a0sph.flatten()])
#     file.write(line + "\n")

# with open('resultados_a0_def.txt', 'w') as a0_file, \
#      open('resultados_da_def.txt', 'w') as da_file:

#     for i in range(N):  # Runge Kutta 4th order

#         phi_ = phi(a0sph)
#         dda_ = dda(a0sph)
#         L1 = h*(da0sph)
#         K1 = h*(dda_)

#         phi_ = phi(a0sph + L1/2)
#         dda_ = dda(a0sph + L1/2)
#         L2 = h*(da0sph + K1/2)
#         K2 = h*(dda_)

#         phi_ = phi(a0sph + L2/2)
#         dda_ = dda(a0sph + L2/2)
#         L3 = h*(da0sph + K2/2)
#         K3 = h*(dda_)

#         phi_ = phi(a0sph + L3)
#         dda_ = dda(a0sph + L3)
#         L4 = h*(da0sph + K3)
#         K4 = h*(dda_)

#         da0sph = da0sph + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
#         a0sph = a0sph + 1/6 * (L1 + 2*L2 + 2*L3 + L4)
        
#         a0_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in a0sph.flatten()])
#         da_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in da0sph.flatten()])
        
#         a0_file.write(a0_line + "\n")
#         da_file.write(da_line + "\n")

r = np.tile(rcol.repeat(PX + 1).reshape(-1, 1), (1, (PR + 1) * (PX + 1)))
x = np.tile(np.tile(xcol, PX + 1).reshape(-1, 1), (1, (PR + 1) * (PX + 1)))

with open('resultados_a_direct5_3.txt', 'w') as a0_file, \
     open('resultados_da_direct5_3.txt', 'w') as da_file:

    for i in range(N):  # Runge Kutta 4th order

        phi = np.dot(a0sph, col.Psi)
        dda = np.dot(np.dot(a0sph, ddrPsi + 2/r*drPsi + (-x**2 + 1)*ddxPsi/r**2 - 2*x*dxPsi), Psi)
        L1 = h*(da0sph)
        K1 = h*(dda)

        phi = np.dot(a0sph + L1/2, col.Psi)
        dda = np.dot(np.dot(a0sph + L1/2, ddrPsi + 2/r*drPsi + (-x**2 + 1)*ddxPsi/r**2 - 2*x*dxPsi), Psi)
        L2 = h*(da0sph + K1/2)
        K2 = h*(dda)

        phi = np.dot(a0sph + L2/2, Psi)
        dda = np.dot(np.dot(a0sph + L2/2,  ddrPsi + 2/r*drPsi + (-x**2 + 1)*ddxPsi/r**2 - 2*x*dxPsi), Psi)
        L3 = h*(da0sph + K2/2)
        K3 = h*(dda)

        phi = np.dot(a0sph + L3, Psi)
        dda = np.dot(np.dot(a0sph + L3,  ddrPsi + 2/r*drPsi + (-x**2 + 1)*ddxPsi/r**2 - 2*x*dxPsi), Psi)
        L4 = h*(da0sph + K3)
        K4 = h*(dda)

        da0sph = da0sph + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
        a0sph = a0sph + 1/6 * (L1 + 2*L2 + 2*L3 + L4)
        
        a0_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in a0sph.flatten()])
        da_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in da0sph.flatten()])
        
        a0_file.write(a0_line + "\n")
        da_file.write(da_line + "\n")


#   t += h
#   # Save the results to the files
#   file_asph.write(str(t))

#   for j in range(len(a0sph)):
#               file_asph.write(str(a0sph[j]))
#               file_asph.write(' ')

#               ##pegar o numero da componente do vetor

#   file_asph.write("\n")

# # Close files
# file_asph.close()