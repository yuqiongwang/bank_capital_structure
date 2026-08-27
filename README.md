Bank Capital Structure Optimization under Basel III
Bayraktar, Chevalier, Ly Vath, Wang
-------------------------------------------------------------

This code solves the bank's optimal dividend, investment, and
recapitalization problem under solvency and liquidity constraints,
formulated as a combined singular and impulse control problem.

-------------------------------------------------------------
Requirements: MATLAB_R2025b

-------------------------------------------------------------
Baseline parameters:

r=0.02, mu=0.04, mu_L=0.03, rho=0.12,
sigma=0.08, sigma_L=0.03, c=0.20, gamma=0.02,
kappa=0.01, kappa_p=0.02, a1=0.045, a2=0.05, a3=0.30.

-------------------------------------------------------------
The scripts should be run in the following order:

1. 00_a1_main_vi_solver.m
Solves the variational inequality and the HJB (without issuance).
Plots the value function v(y), the issuance value \Delta_v,
and the two dimensional value v(x,l). We use the same lower boundary 
data for the HJB case: v^{HJB}(\underline y) = \underline y-1.

2. 00_a2_trajectory_simulation.m
Save the result from a1 first and load it as vi_solution.mat.
This script simulates a single trajectory of Y_t for 50 years.
It gives the plots of accumulative issuance and dividends.

3. 00_b1_sensitivity_analysis.m
Evaluates the value and y^* over grids of (a1, a2, a3).
Reports y* and Delta_v/v(1.2). 

4. 00_b2_constraint_decomposition.m
Compares the value functions under no-Basel, solvency-only,
LCR-only, and joint-constraint regimes, with and without the leverage cap.

5. 00_c1_profit_survival_frontier.m
Solves the regulator's problem: maximize v(y0) subject to
P(\tau >= T) >= \eta, over a full grid of 729 parameter
triples. Computes the Pareto frontier and optimal parameters
for eta = 0.80 and 0.90 (0.85 for $\pi\leq 1$).
OBS: the runtime is large.

6. 00_c2_frontier_with_CI.m
Computes the confidence-adjusted regulatory frontier using Monte
Carlo survival estimates (2000 paths) and 95\% confidence intervals.

7. 00_c3_validate_twosided.m
Independently validate the result and two-sided CI from 00_c2 using 5000 paths.

8. 00_d1_coord_stackelberg_threshold.m
Computes coordinated and bank-preferred intervention thresholds
and the associated recapitalization frequency and issuance outcomes.
OBS: the runtime is very large, around 4 times of c1.

10. 00_d2_kappa_p_sensitivity.m
Studies how the bank's preferred intervention threshold and recapitalization
behavior vary with the proportional issuance cost kappa'.
OBS: the runtime is large.

12. 00_d3_kappa_sensitivity.m
Studies how the bank's preferred intervention threshold and
recapitalization behavior vary with the dilution-cost parameter kappa.
OBS: the runtime is large.




