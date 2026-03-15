%% sensitivity to a_1, a_2, a_3
clear; clc; close all;
%% parameters
base_params.r       = 0.01;
base_params.mu      = 0.04;
base_params.mu_L    = 0.03;
base_params.rho     = 0.12;
base_params.sigma   = 0.08;
base_params.sigma_L = 0.03;
base_params.c       = 0.20;
base_params.gamma   = 0.01;
base_params.a1      = 0.045;   
base_params.a2      = 0.05;
base_params.a3      = 0.30;
base_params.kappa   = 0.01;
base_params.kappa_p = 0.02;
base_params.pi_cap  = [];

settings.y_max       = 2.5;
settings.N           = 301;
settings.dt          = 0.05;
settings.max_iter    = 1500;
settings.tol         = 1e-6;
settings.pi_grid_size = 151;

%evaluation point (need to be in the continuaton)
y_eval = 1.2;  

%% a_1
a1_values = [0.045, 0.05, 0.06, 0.07,0.08, 0.09,0.10, 0.11,0.12];
n_a1 = length(a1_values);
results_a1 = struct();
results_a1.a1 = a1_values;
results_a1.y_star = zeros(n_a1, 1);
results_a1.Delta_v_ratio = zeros(n_a1, 1);

for i = 1:n_a1
    fprintf('  a₁ = %.3f ... ', a1_values(i));
    params = base_params;
    params.a1 = a1_values(i);
    check_params_valid(params);  
    [sol_qvi] = solve_bank_qvi(params, settings, 'QVI');
    [sol_hjb] = solve_bank_qvi(params, settings, 'HJB');
    results_a1.y_star(i) = sol_qvi.y_star;
    v_qvi_eval = interp1(sol_qvi.y, sol_qvi.v, y_eval, 'linear', 'extrap');
    v_hjb_eval = interp1(sol_hjb.y, sol_hjb.v, y_eval, 'linear', 'extrap');
    Delta_v = v_qvi_eval - v_hjb_eval;
    results_a1.Delta_v_ratio(i) = Delta_v / v_qvi_eval;
    end

%% a_2 LCR
%a2_values = [0.05];
a2_values = [0.05,0.06 0.08, 0.09,0.10,0.12, 0.15, 0.18, 0.2];
n_a2 = length(a2_values);
results_a2 = struct();
results_a2.a2 = a2_values;
results_a2.y_star = zeros(n_a2, 1);
results_a2.Delta_v_ratio = zeros(n_a2, 1);
%
for i = 1:n_a2
    fprintf('  a₂ = %.3f ... ', a2_values(i));
    params = base_params;
    params.a2 = a2_values(i);
    check_params_valid(params);
    [sol_qvi] = solve_bank_qvi(params, settings, 'QVI');
    [sol_hjb] = solve_bank_qvi(params, settings, 'HJB');
    
    results_a2.y_star(i) = sol_qvi.y_star;
    
    v_qvi_eval = interp1(sol_qvi.y, sol_qvi.v, y_eval, 'linear', 'extrap');
    v_hjb_eval = interp1(sol_hjb.y, sol_hjb.v, y_eval, 'linear', 'extrap');
    Delta_v = v_qvi_eval - v_hjb_eval;
    results_a2.Delta_v_ratio(i) = Delta_v / v_qvi_eval;
    end

%% a_3 HQLA
%a3_values = [0.15]
a3_values = [0.15, 0.20, 0.25, 0.3,0.35,  0.4 ,0.45, 0.5, 0.55];
n_a3 = length(a3_values);
results_a3 = struct();
results_a3.a3 = a3_values;
results_a3.y_star = zeros(n_a3, 1);
results_a3.Delta_v_ratio = zeros(n_a3, 1);
%
for i = 1:n_a3
    fprintf('  a₃ = %.3f ... ', a3_values(i));
    params = base_params;
    params.a3 = a3_values(i);
    check_params_valid(params);
    [sol_qvi] = solve_bank_qvi(params, settings, 'QVI');
    [sol_hjb] = solve_bank_qvi(params, settings, 'HJB');
    results_a3.y_star(i) = sol_qvi.y_star;
    v_qvi_eval = interp1(sol_qvi.y, sol_qvi.v, y_eval, 'linear', 'extrap');
    v_hjb_eval = interp1(sol_hjb.y, sol_hjb.v, y_eval, 'linear', 'extrap');
    Delta_v = v_qvi_eval - v_hjb_eval;
    results_a3.Delta_v_ratio(i) = Delta_v / v_qvi_eval;
end

%% tables output
%a1
T1 = table(results_a1.a1(:), results_a1.y_star, results_a1.Delta_v_ratio, ...
    'VariableNames', {'a1', 'y_star', 'Delta_v_ratio'});
fprintf('\nSensitivity to a1:\n');
disp(T1);
% a2
T2 = table(results_a2.a2(:), results_a2.y_star, results_a2.Delta_v_ratio, ...
    'VariableNames', {'a2', 'y_star', 'Delta_v_ratio'});
fprintf('\nSensitivity to a2:\n');
disp(T2);
% a3
T3 = table(results_a3.a3(:), results_a3.y_star, results_a3.Delta_v_ratio, ...
    'VariableNames', {'a3', 'y_star', 'Delta_v_ratio'});
fprintf('\nSensitivity to a3:\n');
disp(T3);
