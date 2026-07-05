# Numeric_PDE_illustration
Illustration of numeric methods used to solve one dimensional heat equation and wave equation, their evolution of error over time. 

## Diffusion equation

To solve

$$
u_t = a u_{xx}
$$

$$
0 < x < 1
$$

Given

$$
u(0,t) = u(1,t) = 0,
$$

$$
u(x,0) = \sin(\pi x)
$$

using forward difference in time and central difference for space.

primary code used to apply this is the following

```python
a = 1.0
L = 1.0  # length of space domain
T = 0.1  # Final time

# Step sizes
Nx = 50
dx = L / 50
dt = 0.8 * dx**2 / (2 * a)  # the method is conditionally stable, hence the specific choice for dt
Nt = int(T / dt)

x = np.linspace(0, L, Nx + 1)

r = a * dt / dx**2

# Initial condition
u = np.sin(np.pi * x)

# Boundary condition
u[0] = 0.0
u[-1] = 0.0

for n in range(Nt):
    for i in range(1,Nx):
        u_new[i] =u[i] + r*(u[i-1] - 2*u[i] + u[i+1]) #marching forward in time by applying the scheme
    u[0] = 0.0
    u[-1] = 0.0 # forcing boundry condition

    u[:] = u_new # updating new soltuion
```
