%% Independent Monte Carlo validation 
% of the policies selected by C_2
% c2_frontier_with_CI.m.
% two sided
clear; clc;

%% Baseline parameters
params.r           = 0.02;
params.mu          = 0.04;
params.mu_L        = 0.03;
params.sigma       = 0.08;
params.sigma_L     = 0.03;
params.c           = 0.20;
params.gamma       = 0.02;
params.rho         = 0.12;
params.kappa       = 0.01;
params.kappa_p     = 0.02;
params.underline_y = 1.02;

T_horizon = 5.0;
y0 = 1.20;

settings = struct( ...
    'N', 301, ...
    'y_max', 2.5, ...
    'dt', 0.05, ...
    'max_iter', 1500, ...
    'tol', 1e-6, ...
    'pi_grid_size', 151);

n_paths = 5000;
dt_mc = 0.05;
validation_seed = 20260826;
z_score = 1.96;

%% Policies selected by c2
policies = [
    0.10421, 0.05, 0.25, Inf, 0.80;
    0.12000, 0.05, 0.25, Inf, 0.90;
    0.11605, 0.05, 0.25,   1, 0.80;
    0.12000, 0.05, 0.25,   1, 0.84
];

n = size(policies,1);

%%
regime   = strings(n,1);
eta      = zeros(n,1);
a1       = zeros(n,1);
a2       = zeros(n,1);
a3       = zeros(n,1);
P_val    = zeros(n,1);
P_lo_val = zeros(n,1);
P_hi_val = zeros(n,1);

%%
for i = 1:n
    params.a1 = policies(i,1);
    params.a2 = policies(i,2);
    params.a3 = policies(i,3);
    eta(i)    = policies(i,5);

    if isinf(policies(i,4))
        params.pi_cap = [];
        regime(i) = "pi<inf";
    else
        params.pi_cap = policies(i,4);
        regime(i) = "pi<=1";
    end

    sol = solve_bank_qvi(params, settings, 'QVI');
    rng(validation_seed + i, 'twister');
    P_val(i) = estimate_survival_prob( ...
        sol, params, y0, T_horizon, n_paths, dt_mc);
    x = round(P_val(i) * n_paths);
    [P_lo_val(i), P_hi_val(i)] = wilson_ci( ...
        x, n_paths, z_score);

    a1(i) = params.a1;
    a2(i) = params.a2;
    a3(i) = params.a3;
end

validated = P_lo_val >= eta;

%% Report
T = table( ...
    regime, eta, a1, a2, a3, P_val, P_lo_val, P_hi_val, validated, ...
    'VariableNames', { ...
        'regime','eta','a1','a2','a3', ...
        'P_val','P_lo_val','P_hi_val','validated'});
%% 
function [lo, hi] = wilson_ci(x, n, z)
    phat = x / n;
    denom = 1 + z^2 / n;
    centre = (phat + z^2 / (2*n)) / denom;
    halfw = z * sqrt(phat*(1-phat)/n + z^2/(4*n^2)) / denom;
    lo = max(0, centre - halfw);
    hi = min(1, centre + halfw);
end
