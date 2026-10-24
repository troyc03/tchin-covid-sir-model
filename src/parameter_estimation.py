import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize
from sir_model import sir_model, solve_ivp
import matplotlib.pyplot as plt

# Function to estimate parameters beta and gamma using least squares optimization
def estimate_parameters(t, y_data, initial_conditions):
    def objective(params):
        beta, gamma = params
        sol = solve_ivp(sir_model, (t[0], t[-1]), initial_conditions, args=(beta, gamma), t_eval=t)
        I_model = sol.y[1]  # Infected compartment
        return mean_squared_error(y_data, I_model)

    # Initial guess for beta and gamma
    initial_guess = [0.3, 0.1]
    bounds = [(0, None), (0, None)]  # beta and gamma must be non-negative

    result = minimize(objective, initial_guess, bounds=bounds)
    return result.x  # Returns estimated beta and gamma

if __name__ == "__main__":
    # Example usage
    t = np.linspace(0, 160, 160)  # Time points
    S0 = 0.99
    I0 = 0.01
    R0 = 0.0
    initial_conditions = [S0, I0, R0]

    # Simulated data for testing (replace with real data)
    beta_true = 0.3
    gamma_true = 0.1
    sol = solve_ivp(sir_model, (t[0], t[-1]), initial_conditions, args=(beta_true, gamma_true), t_eval=t)
    y_data = sol.y[1] + np.random.normal(0, 0.01, size=sol.y[1].shape)  # Add noise to the infected data

    estimated_params = estimate_parameters(t, y_data, initial_conditions)
    print(f"Estimated beta: {estimated_params[0]:.4f}, Estimated gamma: {estimated_params[1]:.4f}")

    # Print out equations for the SIR model with estimated parameters
    print(f"SIR Model Equations with Estimated Parameters:")
    print(f"dS/dt = -{estimated_params[0]:.4f} * S * I")
    print(f"dI/dt = {estimated_params[0]:.4f} * S * I - {estimated_params[1]:.4f} * I")
    print(f"dR/dt = {estimated_params[1]:.4f} * I")

    # Plot least squares fit of the SIR model to the data
    sol_estimated = solve_ivp(sir_model, (t[0], t[-1]), initial_conditions, args=(estimated_params[0], estimated_params[1]), t_eval=t)
    plt.plot(t, y_data, 'o', label='Noisy Data', markersize=3)
    plt.plot(t, sol_estimated.y[1], label='Fitted Model', color ='red')
    plt.xlabel('Time (days)')
    plt.ylabel('Proportion of Infected Population')
    plt.title('Least Squares Fit of SIR Model to Data', fontsize=16, fontweight='bold')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()

    # Print R2 score for the fit
    r2_score = 1 - (np.sum((y_data - sol_estimated.y[1])**2) / np.sum((y_data - np.mean(y_data))**2))
    print(f"R2 Score of the Fit: {r2_score:.4f}")