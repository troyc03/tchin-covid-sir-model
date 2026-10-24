import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the SIR model system of differential equations
def sir_model(t, y, beta, gamma):
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

# Plot phase portrait of the SIR model (S vs I plane)
def plot_phase_portrait(beta, gamma, S_range, I_range, num_points=20):
    S = np.linspace(S_range[0], S_range[1], num_points)
    I = np.linspace(I_range[0], I_range[1], num_points)
    S_grid, I_grid = np.meshgrid(S, I)
    
    # Vector field equations
    dSdt = -beta * S_grid * I_grid
    dIdt = beta * S_grid * I_grid - gamma * I_grid
    
    plt.figure(figsize=(12, 8))
    
    # Normalize vectors for cleaner visual layout
    norm = np.hypot(dSdt, dIdt)
    norm[norm == 0] = 1  # Avoid division by zero
    plt.quiver(S_grid, I_grid, dSdt/norm, dIdt/norm, color='blue', alpha=0.3, pivot='mid')
    
    # Plot solution trajectories along the vector field
    for S0 in np.linspace(S_range[0] + 0.1, S_range[1], 5):
        for I0 in np.linspace(I_range[0] + 0.05, I_range[1], 4):
            initial_conditions = [S0, I0, 0.0]
            t_span = (0, 50)
            t_eval = np.linspace(t_span[0], t_span[1], 200)
            sol = solve_ivp(sir_model, t_span, initial_conditions, args=(beta, gamma), t_eval=t_eval)
            plt.plot(sol.y[0], sol.y[1], color='red', alpha=0.6)
            
    # Plot Nullclines
    # 1. dS/dt = 0 nullclines occur at S=0 and I=0
    plt.axhline(0, color='g', linestyle='--', alpha=0.7, label='dS/dt=0 or dI/dt=0 (I=0)')
    
    # 2. dI/dt = 0 nullcline occurs at S = gamma / beta
    S_threshold = gamma / beta
    if S_range[0] <= S_threshold <= S_range[1]:
        plt.axvline(S_threshold, color='k', linestyle='--', alpha=0.7, label=f'dI/dt=0 (S = γ/β = {S_threshold:.2f})')
        
    # Plot critical point (disease-free equilibrium)
    plt.plot(1, 0, 'ko', markersize=8, label='Disease-Free Equilibrium (S=1, I=0)')
    
    plt.xlabel('Susceptible (S)', fontsize=12)
    plt.ylabel('Infected (I)', fontsize=12)
    plt.title(f'Phase Portrait of the SIR Model (beta={beta}, gamma={gamma})', fontsize=16, fontweight='bold')
    plt.xlim(S_range)
    plt.ylim(I_range)
    plt.legend(loc='upper right')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

if __name__ == "__main__":
    # Parameters
    beta = 0.3
    gamma = 0.1
    
    # Compute and print eigenvalues
    eigenvalues = compute_eigenvalues(beta, gamma)
    print(f"Eigenvalues at disease-free equilibrium: {eigenvalues}")
    
    # Plot phase portrait
    S_range = (0, 1.1)
    I_range = (0, 0.6)
    plot_phase_portrait(beta, gamma, S_range, I_range)
