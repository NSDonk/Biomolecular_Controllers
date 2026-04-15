# Biomolecular_Controllers
Code for modeling proportional, proportional-integral, proportional-derivative and proportional-integral-derivative controllers. 

## Overview
This project focuses on the modeling and analysis of biomolecular feedback control systems as dynamical systems through a control theory analytical view. The goal is to understand how classical control architectures (proportional, proportional-integral, proportional-derivative, and proportional-integral-derivative) can be implemented in biochemical reaction networks and used to regulate biological processes.

---

## Methods
- **Dynamical Systems Modeling**
  - Ordinary differential equation (ODE) models of biochemical reaction networks
  - Nonlinear system simulation and steady-state analysis

- **Control-Theoretic Analysis**
  - Mapping biomolecular circuits to classical controller architectures (P, PI, PD, PID)
  - Derivation and analysis of transfer functions
  - Evaluation of steady-state error, stability, and transient dynamics

- **Simulation and Parameter Exploration**
  - Parameter sweeps to evaluate controller performance across regimes
  - Analysis of stability conditions and sensitivity to gain parameters

---

## Results
- Demonstrated how biomolecular implementations of classical controllers exhibit characteristic tradeoffs between stability, responsiveness, and steady-state error 
- Identified conditions under which integral action enables robust setpoint tracking, as well as regimes where system dynamics lead to degraded performance  
- Explored how nonlinear biochemical constraints shape the achievable behavior of feedback controllers  

---

## Notes
This project is part of ongoing work in biomolecular control systems and synthetic biology. Future work includes extending the framework to other complex biological systems as this notebook currently focuses on modeling thyroid disease (Hashimoto's and Graves' disease, hypo- / hyperthyroidism respectively).
