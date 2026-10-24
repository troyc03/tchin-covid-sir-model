# Data-Driven Modeling of COVID-19 Dynamics Using the SIR Model

## Overview

This project investigates the spread of COVID-19 through data-driven mathematical modeling. Using historical COVID-19 case data, the project applies the Susceptible-Infectious-Removed (SIR) model to study epidemic dynamics, estimate model parameters, and analyze the evolution of infections over time.

The main objective is to determine how effectively a system of differential equations can describe the observed progression of a real-world epidemic. Both the classical nonlinear SIR model and its linearized approximation will be explored to understand early epidemic growth, infection peaks, and long-term behavior.

## Research Questions

- How accurately can the SIR model describe the observed spread of COVID-19?
- What transmission and removal rates can be estimated from reported case data?
- How does the linearized SIR model approximate the early exponential growth of an epidemic?
- When do the susceptible, infectious, and removed populations approach a steady state?
- How sensitive are the model's predictions to changes in its parameters?
- What are the limitations of using a simplified SIR model to represent COVID-19 dynamics?

## Mathematical Framework

The classical SIR model divides the population into three compartments:

- **Susceptible (S):** Individuals who can become infected.
- **Infectious (I):** Individuals who can transmit the disease.
- **Removed (R):** Individuals who no longer contribute to transmission.

The model is governed by the following system of ordinary differential equations:

$\dot{S} = \frac{dS}{dt}=-\beta\frac{SI}{N}$

$\dot{I} = \frac{dI}{dt}=\beta\frac{SI}{N}-\gamma I$

$\dot{R} = \frac{dR}{dt}=\gamma I$

where \(N\) is the total population, \(\beta\) is the transmission rate, and \(\gamma\) is the removal rate.

The basic reproduction number is defined as

$\mathcal{R}_0=\frac{\beta}{\gamma}$

A linearized approximation will also be investigated to characterize early epidemic growth and estimate the exponential growth rate and doubling time.

## Methodology

The project will follow these steps:

1. **Exploratory Data Analysis:** Examine cumulative COVID-19 cases, calculate changes in reported cases, identify epidemic waves, and investigate reporting irregularities.
2. **Data Preprocessing:** Select a geographic region and study period, account for the dataset's reporting cadence, and prepare observations for model fitting.
3. **Model Development:** Formulate the nonlinear SIR model and derive its linearized approximation for early epidemic dynamics.
4. **Parameter Estimation:** Estimate model parameters using numerical optimization and appropriate fitting methods.
5. **Numerical Simulation:** Solve the differential equations numerically and simulate the evolution of the three population compartments.
6. **Model Validation:** Compare model predictions with observed case data and evaluate prediction errors.
7. **Dynamical Analysis:** Investigate epidemic growth, the infectious peak, long-term equilibrium, and sensitivity to parameter changes.
8. **Limitations and Interpretation:** Assess the effects of simplified assumptions, incomplete case reporting, and changes in transmission conditions.

## Dataset

The project uses historical COVID-19 cumulative confirmed-case data covering multiple geographic regions over several years.

Because cumulative confirmed cases do not directly measure the number of currently infectious individuals, an observation model will be needed to connect reported cases to the SIR compartments. Reporting delays, underreporting, and irregular updates will also be considered during data preprocessing and interpretation.

## Tools and Technologies

- **Python** for data processing, analysis, and simulation
- **Pandas and NumPy** for data manipulation and numerical calculations
- **Matplotlib** for visualizing epidemic trends and model predictions
- **SciPy** for numerical integration and parameter estimation
- **Differential equations and dynamical systems theory** for mathematical analysis

## Expected Outcomes

This project aims to:

- Develop a mathematical model of COVID-19 transmission from observed data.
- Estimate epidemiological parameters and investigate their interpretation.
- Compare linearized and nonlinear model behavior.
- Evaluate how well the SIR model reproduces observed epidemic trends.
- Estimate when the modeled epidemic approaches a practical steady state.
- Identify the limitations of the model and potential improvements.

The results will demonstrate how mathematical modeling, numerical simulation, and real-world data can be combined to investigate the dynamics of infectious diseases.

## Project Status

**In progress.** Current work focuses on exploratory data analysis, deriving the SIR model, developing an analytical linearized approximation, and preparing the data for parameter estimation and numerical simulation.
