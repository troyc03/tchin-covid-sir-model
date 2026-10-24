import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from sir_model import sir_model

# Perform stability analysis by calculating the Jacobian matrix at the disease-free equilibrium
def jacobian_matrix(beta, gamma):
    # Disease-free equilibrium: S=1, I=0, R=0
    J = np.array([[-beta, 0, 0],
                  [beta, -gamma, 0],
                  [0, gamma, 0]])
    return J

# Function to compute eigenvalues of the Jacobian matrix
def compute_eigenvalues(beta, gamma):
    J = jacobian_matrix(beta, gamma)
    eigenvalues = np.linalg.eigvals(J)
    return eigenvalues

# Plot phase portrait of the SIR model
def plot_phase_portrait(beta, gamma, S_range, I_range, num_points=20):
    S = np.linspace(S_range[0], S_range[1], num_points)
    I = np.linspace(I_range[0], I_range[1], num_points)
    S_grid, I_grid = np.meshgrid(S, I)

    dSdt = -beta * S_grid * I_grid
    dIdt = beta * S_grid * I_grid - gamma * I_grid

    plt.figure(figsize=(12, 8))
    plt.quiver(S_grid, I_grid, dSdt, dIdt, color='blue', alpha=0.5) 

    # Plot solution trajectories along the vector field
    for S0 in np.linspace(S_range[0], S_range[1], 5):
        for I0 in np.linspace(I_range[0], I_range[1], 5):
            initial_conditions = [S0, I0, 0.0]
            t_span = (0, 50)
            t_eval = np.linspace(t_span[0], t_span[1], 200)
            sol = solve_ivp(sir_model, t_span, initial_conditions, args=(beta, gamma), t_eval=t_eval)
            plt.plot(sol.y[0], sol.y[1], color='red')  # Plot S vs I

    # Plot critical point (disease-free equilibrium)
    plt.plot(1, 0, 'ko', label='Disease-Free Equilibrium (S=1, I=0)')

    plt.ylabel('Infected (I)')
    plt.title(f'Phase Portrait of the SIR Model (beta={beta}, gamma={gamma})', fontsize=16, fontweight='bold')
    plt.xlim(S_range)
    plt.ylim(I_range)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

if __name__ == "__main__":
    # Example usage
    beta = 0.3
    gamma = 0.1

    # Compute and print eigenvalues
    eigenvalues = compute_eigenvalues(beta, gamma)
    print(f"Eigenvalues of the Jacobian matrix at disease-free equilibrium: {eigenvalues}")

    # Plot phase portrait
    S_range = (0, 1)
    I_range = (0, 0.5)
    plot_phase_portrait(beta, gamma, S_range, I_range)
