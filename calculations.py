import numpy as np
from variables import Nx, Nt, dt, dx, theta_new, theta_old, F, A, B, u_array

# Copy initial state for the 'new' array for plotting
theta_new[:] = theta_old[:]
theta_initial = np.copy(theta_old)

# Print to keep user informed
print("Starting implicit advection simulation")
print(f"Total time time steps: {Nt}")

# Check CFL condition
CFL = np.max(np.abs(u_array) * dt/dx)
print(f"CFL number: {CFL:.3f}")

current_time = 0.0

# Time iteration loop (j=1 to Nt)
for j in range(1, Nt + 1):
    current_time = j * dt

    # Assume fixed boundary conditions so:
    theta_new[0] = theta_initial[0]

    # Calculate right hand side vector F
    F[:] = theta_old[1:]

    # Solve the linear system for I = 1 to Nx
    for I in range(1, Nx):
        k = I - 1 
        theta_new[I] = (F[k] - B[k] * theta_new[I - 1]) / A[k]

theta_old[:] = theta_new[:]

print(f"Simulation complete at t = {current_time:.1f} s")