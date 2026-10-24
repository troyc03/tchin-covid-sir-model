import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the SIR model system of differential equations
def linearized_sir_model(t, y, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

# Calculate the Jacobian matrix at the disease-free equilibrium (S=1, I=0, R=0)
def jacobian_matrix(beta, gamma):
    # Formed by partial derivatives evaluated at S=1, I=0, R=0
    J = np.array([
        [0, -beta, 0],
        [0, beta - gamma, 0],
        [0, gamma, 0]
    ])
    return J

# Compute eigenvalues of the Jacobian matrix
def compute_eigenvalues(beta, gamma):
    J = jacobian_matrix(beta, gamma)
    eigenvalues = np.linalg.eigvals(J)
    return eigenvalues

# Run baseline simulation
sol = solve_ivp(linearized_sir_model, [0, 50], [0.99, 0.01, 0.0], args=(0.3, 0.1), t_eval=np.linspace(0, 50, 200))
eigenvalues = compute_eigenvalues(0.3, 0.1)

# Plot dynamics over time
plt.figure(figsize=(12, 6))
plt.plot(sol.t, sol.y[0], label='Susceptible (S)', color='blue')
plt.plot(sol.t, sol.y[1], label='Infected (I)', color='red')
plt.plot(sol.t, sol.y[2], label='Recovered (R)', color='green')
plt.title('Linearized SIR Model Dynamics', fontsize=16, fontweight='bold')
plt.xlabel('Time')
plt.ylabel('Population Fraction')
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.show()
