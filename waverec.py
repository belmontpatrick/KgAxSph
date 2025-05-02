import numpy as np
from initialdata import gaussian_approx
from parameters import N, px, A0
from initialdata import Psi_plot,Rplot,Xplot
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
# import scipy.special as sp

data = np.loadtxt('resultados_a_def_tst.txt')
# data01 = data[99,1:]
# print(data01)
 
def gaussian_approx(r,x,b0):
    res = sum(b0[k] * Psi_plot(k,r,x) for k in range((N+1) * (px+1)))
    return res

# phi = gaussian_approx(Rplot,Xplot,data01)

# fig3 = plt.figure()
# ax3 = plt.axes(projection = '3d')
# ax3.plot_surface(Rplot, Xplot, phi, cmap='viridis')
# ax3.set_xlabel('r')
# ax3.set_ylabel('x')
# ax3.set_zlabel('phi')
# # plt.savefig('exact_phi_trunc121.png', dpi=300, bbox_inches='tight')
# plt.show()


###Animation

# plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'  #Ajuste o caminho se necessário
# get_ipython().run_line_magic('matplotlib', 'notebook')

# data = np.loadtxt('resultados_a0.txt')
times = data[:, 0]  #time column
all_a0 = data[:, 1:]  #a0 
r_min, r_max = 0.1, 5.0
x_min, x_max = -1, 1


#grid
r_vals = np.linspace(r_min, r_max, 30)
x_vals = np.linspace(x_min, x_max, 30)
R, X = np.meshgrid(r_vals, x_vals)

#psi para plot 3d
total_coeffs = (N+1)*(px+1)
Psi_grid = np.zeros((total_coeffs, *R.shape))

for k in range(total_coeffs):
    # Vectorize manualmente para cada k
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            Psi_grid[k,i,j] = Psi_plot(k, R[i,j], X[i,j])

#phi
phi_evolution = np.zeros((len(times), *R.shape))

for t_idx in range(len(times)):
    #Produto vetorizado: a0[t_idx] • Psi_grid
    phi_evolution[t_idx] = np.sum(all_a0[t_idx][:,None,None] * Psi_grid, axis=0)


#3D animation
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('r')
ax.set_ylabel('x')
ax.set_zlabel('ϕ(r,x,t)')
ax.set_zlim(np.min(phi_evolution), np.max(phi_evolution))

#initial frame
t_idx = 0
surf = ax.plot_surface(R, X, phi_evolution[t_idx], 
                      cmap='viridis', 
                      rstride=1, cstride=1)
ax.set_title(f'Tempo = {times[t_idx]:.3f}')

#atualização da função
def update(frame):
    global surf
    
    t_idx = min(frame*5, len(times)-1)
    
    if surf:
        surf.remove()
    
    surf = ax.plot_surface(R, X, phi_evolution[t_idx], 
                         cmap='viridis',
                         rstride=1, cstride=1)
    ax.set_title(f'Tempo = {times[t_idx]:.3f}')
    return surf

#animação
ani = FuncAnimation(fig, update,
                    frames=len(times)//5,
                    # frames=len(times),
                    interval=200,  
                    blit=False,   
                    repeat=True)

plt.tight_layout()
plt.show()

"#save it as a gif"
# ani.save("evolucao_phi_def.gif", writer='pillow', fps=15, dpi=100)

"save it as mp4"
ani.save("evolucao_phi_def.mp4", writer='ffmpeg',  # Changed from 'pillow' to 'ffmpeg'
         fps=13, 
         dpi=100,
         bitrate=1800)