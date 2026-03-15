%% Optimized bank
clear; clc; close all;
%% parameters
params = struct();
params.r       = 0.01;
params.mu      = 0.04;
params.mu_L    = 0.03;
params.sigma   = 0.08;
params.sigma_L = 0.03;
params.c       = 0.20;
params.gamma   = 0.01;
params.rho     = 0.12;
params.kappa   = 0.01;
params.kappa_p = 0.02;
params.pi_cap  = 1;% either give [] or 1 here

params.a1 = 0.12;% use the optimized one from a4 
params.a2 = 0.05;
params.a3 = 0.25;

settings = struct();
settings.N = 301;
settings.y_max = 4;
settings.dt = 0.05;
settings.max_iter = 1500;
settings.tol = 1e-6;
settings.pi_grid_size = 151;

n_paths = 1000;
T_sim = 50;
dt_sim = 0.05;

y0_values = [1.05, 1.10, 1.15, 1.20, 1.25, 1.30];
n_banks = length(y0_values);
%% solve the bank variational inequality
sol = solve_bank_qvi(params, settings, 'QVI');
y_star = sol.y_star;

y = sol.y;
v = sol.v;
v_prime = sol.v_prime;
pi_star_grid = sol.pi_star;
[y_post, y_post_method] = pick_y_post(y, v, v_prime, y_star, params.kappa_p);
dy = y(2) - y(1);
eps_hit = max(1e-6, 0.5*dy);
y_trigger = 1 + eps_hit;
%% MC simulation
n_steps = ceil(T_sim / dt_sim);
sqrt_dt = sqrt(dt_sim);
rng(42);
dW_all = sqrt_dt * randn(n_paths, n_steps);
%
results = struct();
results.y0 = y0_values;
results.y_star = y_star;
results.y_post = y_post;
results.E_total_xi = zeros(n_banks, 1);    
results.E_total_dZ = zeros(n_banks, 1);    
results.Sharpe = zeros(n_banks, 1);
for bank = 1:n_banks
    y0 = y0_values(bank);
    total_xi_path = zeros(n_paths, 1);
    total_dZ_path = zeros(n_paths, 1);
    iss_cost_path = zeros(n_paths, 1);
    for path = 1:n_paths
        Y = y0;
        total_xi = 0;
        total_dZ = 0;
        total_cost = 0;
        for step = 1:n_steps
            pi_opt = interp1(y, pi_star_grid, Y, 'linear', 'extrap');
            pi_opt = max(0, min(compute_pi_bar(Y, params), pi_opt));
            [b, sigma2] = drift_and_sigma2(Y, pi_opt, params);
            sigma_Y = sqrt(max(0, sigma2));
            dW = dW_all(path, step);
            Y_new = Y + b * dt_sim + sigma_Y * dW;
            if Y_new >= y_star
                dZ = Y_new - y_star;
                total_dZ = total_dZ + dZ;
                Y_new = y_star;
            end
            if Y_new <= y_trigger
                y_before = Y_new;
                xi = (y_post - y_before) / (1 - params.kappa_p);
                total_xi = total_xi + xi;
                total_cost = total_cost + params.kappa + params.kappa_p * xi;
                Y_new = y_post;
            end
            Y = max(1.0, Y_new);
        end
        
        total_xi_path(path) = total_xi;
        total_dZ_path(path) = total_dZ;
        iss_cost_path(path) = total_cost;
    end
    net_payout = total_dZ_path - iss_cost_path;
    results.E_total_xi(bank) = mean(total_xi_path);
    results.E_total_dZ(bank) = mean(total_dZ_path);
    std_net = std(net_payout);
    %sharpe ratio
    if std_net > 1e-10
        results.Sharpe(bank) = mean(net_payout) / std_net;
    else
        results.Sharpe(bank) = NaN;
    end
end

%% table
for bank = 1:n_banks
    fprintf('%.2f   | %10.4f | %10.4f | %10.4f\n', ...
        results.y0(bank), ...
        results.E_total_xi(bank), ...
        results.E_total_dZ(bank), ...
        results.Sharpe(bank));
end