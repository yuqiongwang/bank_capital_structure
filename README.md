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

r=0.01, mu=0.04, mu_L=0.03, rho=0.12,
sigma=0.08, sigma_L=0.03, c=0.20, gamma=0.01,
kappa=0.01, kappa_p=0.02, a1=0.045, a2=0.05, a3=0.30.

-------------------------------------------------------------
The scripts should be run in the following order:

1. a1_main_vi_solver.m
Solves the variational inequality and the HJB (without issuance).
Plots the value function v(y), the issuance value \Delta_v,
and the two dimensional value v(x,l). This is Figure 1 in the paper.

2. a2_trajectory_simulation.m
Save the result from a1 first and load it as vi_solution.mat.
This script simulates a single trajectory of Y_t for 50 years.
It gives the plots of accumulative issuance and dividends.
They correspond to Figure 2 and 3 in the paper.

3. a3_sensitivity_analysis.m
Evaluates the value and y^* over grids of (a1, a2, a3).
Reports y* and Delta_v/v(1.2). This is Table 1 in the paper.

4. a4_health_constrained_opt.m
Solves the regulator's problem: maximize v(y0) subject to
P(\tau >= T) >= \eta, over a full grid of 729 parameter
triples. Computes the Pareto frontier and optimal parameters
for eta = 0.80 and 0.90. It yields Tables 2-5 and Figure 4
in the paper. OBS: the runtime is large.

5. a5_optimized_trajectory.m
Uses the optimized parameters from a4 to simulate 1000 Monte
Carlo paths for 6 banks with different starting points.
Reports expected cumulative issuance, dividends, and the Sharpe
ratio of the net gain. It yields Tables 6-7 of the paper.
