import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy.special as sp

# Load data
data = np.loadtxt('resultados_a0.txt')
times = data[:, 0]  # time column
all_a0 = data[:, 1:]  # a0 coefficients

# Define your parameters (set these to your actual values)
N = 10  # Example value
px = 10  # Example value

# Define your Psi_plot function
def Psi_plot(k, r, x):
    """Your actual Psi_plot function implementation"""
    return sp.jv(0, r) * np.cos(x)  # Example placeholder

# Grid parameters
r_min, r_max = 0.1, 5.0
x_min, x_max = -1, 1

# Create grid
r_vals = np.linspace(r_min, r_max, 30)
x_vals = np.linspace(x_min, x_max, 30)
R, X = np.meshgrid(r_vals, x_vals)

# Pre-compute Psi grid
total_coeffs = (N+1)*(px+1)
Psi_grid = np.zeros((total_coeffs, *R.shape))

for k in range(total_coeffs):
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            Psi_grid[k,i,j] = Psi_plot(k, R[i,j], X[i,j])

# Compute phi evolution
phi_evolution = np.zeros((len(times), *R.shape))
for t_idx in range(len(times)):
    phi_evolution[t_idx] = np.sum(all_a0[t_idx][:,None,None] * Psi_grid, axis=0)

# Create contourf animation
fig, ax = plt.subplots(figsize=(12, 8))

# Set up initial contour
contour = ax.contourf(R, X, phi_evolution[0], levels=20, cmap='viridis')
plt.colorbar(contour, label='ϕ(r,x,t)')
ax.set_xlabel('r')
ax.set_ylabel('x')
ax.set_title(f'Tempo = {times[0]:.3f}')

# Update function for animation
def update(frame):
    global contour
    
    t_idx = min(frame, len(times)-1)  # Remove the *5 if you want every timestep
    
    # Clear previous contours
    for c in contour.collections:
        c.remove()
    
    # Create new contour
    contour = ax.contourf(R, X, phi_evolution[t_idx], levels=20, cmap='viridis')
    ax.set_title(f'Tempo = {times[t_idx]:.3f}')
    return contour.collections

# Create animation
ani = FuncAnimation(fig, update,
                   frames=len(times),  # Using all frames now
                   interval=200,  
                   blit=True,   # Changed to True for contourf
                   repeat=True)

plt.tight_layout()
plt.show()

# Save as GIF
ani.save("evolucao_phi_contourf.gif", writer='pillow', fps=15, dpi=100)

# To save as MP4 (requires ffmpeg)
# ani.save("evolucao_phi_contourf.mp4", writer='ffmpeg', fps=15, dpi=100)