import numpy as np
import basis as bs
import initial_data as ini
import timegrid as gr
import dynsys as ds
import parameters as par

np.set_printoptions(precision=16)

PR = par.PR
PX = par.PX
N = gr.N
a0sph = ini.a0sph
da0sph = ini.da0sph
h = gr.h
#t = gr.t
N = gr.N
phi = ds.phi
dda = ds.dda

#Runge-Kutta 4th order


file_asph = open("rk4_asph.txt", "w")
# file_dacyl = open("rk4_dacyl.txt", "w")

t = 0.0           # Initial time

for i in range(N):  # Runge Kutta 4th order

  phi_ = phi(a0sph)
  dda_ = dda(a0sph)
  L1 = h*da0sph
  K1 = h*(dda_)

  phi_ = phi(a0sph + L1/2)
  dda_ = dda(a0sph + L1/2)
  L2 = h*(da0sph + K1/2)
  K2 = h*(dda_)

  phi_ = phi(a0sph + L2/2)
  dda_ = dda(a0sph + L2/2)
  L3 = h*(da0sph + K2/2)
  K3 = h*(dda_)

  phi_ = phi(a0sph + L3)
  dda_ = dda(a0sph + L3)
  L4 = h*(da0sph + K3)
  K4 = h*(dda_)

  da0sph = da0sph + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
  a0sph = a0sph + 1/6 * (L1 + 2*L2 + 2*L3 + L4)

  t += h
  # Save the results to the files
  file_asph.write(str(t))

  for j in range(len(a0sph)):
              file_asph.write(str(a0sph[j]))
              file_asph.write(' ')

  file_asph.write("\n")

# Close files
file_asph.close()