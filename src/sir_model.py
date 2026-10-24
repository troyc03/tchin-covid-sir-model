import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

def sir_model(t, y, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]


if __name__ == "__main__":
    # Example usage
    beta = 0.3
    gamma = 0.1
    S0 = 0.99
    I0 = 0.01
    R0 = 0.0
    t_span = (0, 160)
    t_eval = np.linspace(t_span[0], t_span[1], 160)  # 1-dimensional array of time points for evaluation

    initial_conditions = [S0, I0, R0]
    sol = solve_ivp(sir_model, t_span, initial_conditions, args=(beta, gamma), t_eval=t_eval)

    plt.plot(sol.t, sol.y[0], label='Susceptible', color='black')
    plt.plot(sol.t, sol.y[1], label='Infected', color='red')
    plt.plot(sol.t, sol.y[2], label='Recovered', color='green')
    plt.xlabel('Time (days)')
    plt.ylabel('Proportion of Population')
    plt.title('SIR Model Simulation', fontsize=16, fontweight='bold')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()