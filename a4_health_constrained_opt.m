%% constrained optimization problem for the regulator
clear; clc; close all;
%% parameters (model)
base_params = struct();
base_params.r = 0.01;           
base_params.mu = 0.04;          
base_params.sigma = 0.08;       
base_params.mu_L = 0.03;       
base_params.sigma_L = 0.03;     
base_params.c = 0.20;           
base_params.gamma = 0.01;       
base_params.rho = 0.12;        
base_params.kappa = 0.01;       
base_params.kappa_p = 0.02;     
base_params.pi_cap = 1;       

% regulatory parameters initial
base_params.a1 = 0.045;         
base_params.a2 = 0.05;          
base_params.a3 = 0.30;          
% keep the grid same as the sensitivity, to be consistent
a1_grid = [0.045, 0.05, 0.06, 0.07,0.08, 0.09,0.10, 0.11,0.12];   
a2_grid = [0.05,0.06 0.08, 0.09,0.10,0.12, 0.15, 0.18, 0.2];     
a3_grid = [0.15, 0.20, 0.25, 0.3,0.35,  0.4 ,0.45, 0.5, 0.55];    
eta_values = [0.80, 0.9];  % probability threshold
T_horizon = 5.0;  
y0 = 1.20; 
% the following setting should be the same as the main solver
% otherwise the results would differ
settings = struct();
settings.N = 301;
settings.y_max = 2.5;   
settings.dt = 0.05; 
settings.max_iter = 1500;
settings.tol = 1e-6;
settings.pi_grid_size = 151;  
%check parameters
check_params_valid(base_params);
%simulation parameters
n_paths = 1000;% MC number for the probability
dt = 0.05;
%% grid
n_a1 = length(a1_grid);
n_a2 = length(a2_grid);
n_a3 = length(a3_grid);
n_total = n_a1 * n_a2 * n_a3;
%store results
results = struct();
results.a1 = zeros(n_total, 1);
results.a2 = zeros(n_total, 1);
results.a3 = zeros(n_total, 1);
results.y_star = zeros(n_total, 1);
results.v_y0 = zeros(n_total, 1);         
results.P_survival = zeros(n_total, 1);     % P(\tau \geq T)

idx = 0;
for i1 = 1:n_a1
    for i2 = 1:n_a2
        for i3 = 1:n_a3
            idx = idx + 1;
            
            params = base_params;
            params.a1 = a1_grid(i1);
            params.a2 = a2_grid(i2);
            params.a3 = a3_grid(i3);
            
            try
                [sol] = solve_bank_qvi(params, settings, 'QVI');
                v_y0 = interp1(sol.y, sol.v, y0, 'linear', 'extrap');
                P_surv = estimate_survival_prob(sol, params, y0, T_horizon, n_paths, dt);           
                results.a1(idx) = params.a1;
                results.a2(idx) = params.a2;
                results.a3(idx) = params.a3;
                results.y_star(idx) = sol.y_star;
                results.v_y0(idx) = v_y0;
                results.P_survival(idx) = P_surv;    
            catch ME
                results.a1(idx) = params.a1;
                results.a2(idx) = params.a2;
                results.a3(idx) = params.a3;
                results.y_star(idx) = NaN;
                results.v_y0(idx) = NaN;
                results.P_survival(idx) = NaN;
            end
        end
    end
end

%% solve optimization
opt_results = struct();
opt_results.eta = eta_values;
opt_results.a1_opt = zeros(length(eta_values), 1);
opt_results.a2_opt = zeros(length(eta_values), 1);
opt_results.a3_opt = zeros(length(eta_values), 1);
opt_results.v_opt = zeros(length(eta_values), 1);
opt_results.P_opt = zeros(length(eta_values), 1);
opt_results.y_star_opt = zeros(length(eta_values), 1);
opt_results.feasible_count = zeros(length(eta_values), 1);

for k = 1:length(eta_values)
    eta = eta_values(k);
    feasible = results.P_survival >= eta & ~isnan(results.v_y0);
    n_feasible = sum(feasible);
    opt_results.feasible_count(k) = n_feasible;
    
    if n_feasible > 0%find then maximum v0
        v_feasible = results.v_y0;
        v_feasible(~feasible) = -Inf;
        [v_max, idx_max] = max(v_feasible);
        
        opt_results.a1_opt(k) = results.a1(idx_max);
        opt_results.a2_opt(k) = results.a2(idx_max);
        opt_results.a3_opt(k) = results.a3(idx_max);
        opt_results.v_opt(k) = results.v_y0(idx_max);
        opt_results.P_opt(k) = results.P_survival(idx_max);
        opt_results.y_star_opt(k) = results.y_star(idx_max);
    else
        opt_results.a1_opt(k) = NaN;
        opt_results.a2_opt(k) = NaN;
        opt_results.a3_opt(k) = NaN;
        opt_results.v_opt(k) = NaN;
        opt_results.P_opt(k) = NaN;
        opt_results.y_star_opt(k) = NaN;
    end
end

%% compute the pareto frontier(maximizing both)
valid = ~isnan(results.v_y0) & ~isnan(results.P_survival);
v_valid = results.v_y0(valid);
P_valid = results.P_survival(valid);
idx_valid = find(valid);
is_pareto = false(length(v_valid), 1);
for i = 1:length(v_valid)
    dominated = false;% there is no point dominating it, thus pareto
    for j = 1:length(v_valid)
        if i ~= j
            if v_valid(j) >= v_valid(i) && P_valid(j) >= P_valid(i) && ...
               (v_valid(j) > v_valid(i) || P_valid(j) > P_valid(i))
                dominated = true;
                break;
            end
        end
    end
    is_pareto(i) = ~dominated;
end
pareto_idx = idx_valid(is_pareto);
n_pareto = length(pareto_idx);
% sort the list
[~, sort_idx] = sort(results.P_survival(pareto_idx));
pareto_sorted = pareto_idx(sort_idx);
for i = 1:n_pareto
    idx = pareto_sorted(i);
    fprintf('%.3f  %.2f   %.2f   | %.3f    %.4f   %.3f\n', ...
        results.a1(idx), results.a2(idx), results.a3(idx), ...
        results.y_star(idx), results.v_y0(idx), results.P_survival(idx));
end
%% tables
T_pareto = table(...
    results.a1(pareto_sorted), ...
    results.a2(pareto_sorted), ...
    results.a3(pareto_sorted), ...
    results.y_star(pareto_sorted), ...
    results.v_y0(pareto_sorted), ...
    results.P_survival(pareto_sorted), ...
    'VariableNames', {'a1','a2','a3','y_star','v_y0','P_survival'});
fprintf('\nPareto frontier (sorted by survival probability):\n');
disp(T_pareto);
% th optimal one
T_opt = table(...
    eta_values(:), ...
    opt_results.a1_opt, ...
    opt_results.a2_opt, ...
    opt_results.a3_opt, ...
    opt_results.y_star_opt, ...
    opt_results.v_opt, ...
    opt_results.P_opt, ...
    'VariableNames', {'eta','a1_opt','a2_opt','a3_opt','y_star','v_y0','P_survival'});
fprintf('\nOptimal regulatory parameters by safety threshold:\n');
disp(T_opt);

%% plot the frontier (scatter)
figure(1)
scatter(results.P_survival(valid), results.v_y0(valid), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
scatter(results.P_survival(pareto_idx), results.v_y0(pareto_idx), 80, 'r', 'filled', 'MarkerEdgeColor', 'k');
xlabel('P(\tau \geq T)');
ylabel('v(y_0)');
title('pareto');
legend('All configs', 'Pareto frontier', 'Location', 'best');
grid on;
for k = 1:length(eta_values)
    xline(eta_values(k), '--', sprintf('\\eta=%.2f', eta_values(k)), 'LabelVerticalAlignment', 'bottom');
end
