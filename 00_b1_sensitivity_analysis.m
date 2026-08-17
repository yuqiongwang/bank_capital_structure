%% sensitivity to a_1, a_2, a_3
clear; clc; close all;

%% parameters
base_params.r       = 0.02;
base_params.mu      = 0.04;
base_params.mu_L    = 0.03;
base_params.rho     = 0.12;
base_params.sigma   = 0.08;
base_params.sigma_L = 0.03;
base_params.c       = 0.20;
base_params.gamma   = 0.02;

base_params.a1      = 0.045;
base_params.a2      = 0.05;
base_params.a3      = 0.30;

base_params.kappa   = 0.01;
base_params.kappa_p = 0.02;
base_params.pi_cap  = [];        % [] unrestricted; use 1 for leverage cap
base_params.underline_y = 1.02;

settings.y_max        = 2.5;
settings.N            = 301;
settings.dt           = 0.05;
settings.max_iter     = 1500;
settings.tol          = 1e-6;
settings.pi_grid_size = 151;
settings.verbose      = false;
settings.compute_analytical_ystar = false;

% evaluation point
y_eval = 1.20;

%% diagnostics
r_L = base_params.mu_L - base_params.gamma;
if base_params.r <= r_L
    warning('r > r_L violated: r = %.4g, r_L = %.4g.', base_params.r, r_L);
else
    fprintf('condition r > r_L: ok (r=%.4f, r_L=%.4f)\n', base_params.r, r_L);
end

a3_crit = (base_params.mu - base_params.r) / (base_params.rho - base_params.r);
fprintf('Unrestricted feasible floor for a3 when a3>a1: a3 > %.4f\n', a3_crit);

[ok0, info0] = assumption31_holds(base_params);
if ~ok0
    error('violates Assumption 3.1: rho=%.4f, threshold=%.4f.', ...
        base_params.rho, info0.threshold);
end

%% sensitivity to a1
a1_values = [0.045, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12];

results_a1 = table();
results_a1.a1 = a1_values(:);
results_a1.y_star = NaN(numel(a1_values),1);
results_a1.Delta_v_ratio = NaN(numel(a1_values),1);

for i = 1:numel(a1_values)
    p = base_params;              % IMPORTANT: reset baseline
    p.a1 = a1_values(i);

    fprintf('  a1 = %.3f ... ', p.a1);

    [ok, info] = assumption31_holds(p);
    if ~ok
        fprintf('infeasible, check');
        continue;
    end

    check_params_valid(p);
    [y_star, v_qvi_eval, ratio] = compute_one_sensitivity(p, settings, y_eval);

    results_a1.y_star(i) = y_star;
    results_a1.v_qvi_eval(i) = v_qvi_eval;
    results_a1.Delta_v_ratio(i) = ratio;

    fprintf('done: y*=%.4f, ratio=%.4f\n', y_star, ratio);
end

%% sensitivity to a2
a2_values = [0.05, 0.06, 0.08, 0.09, 0.10, 0.12, 0.15, 0.18, 0.20];

results_a2 = table();
results_a2.a2 = a2_values(:);
results_a2.y_star = NaN(numel(a2_values),1);
results_a2.Delta_v_ratio = NaN(numel(a2_values),1);

for i = 1:numel(a2_values)
    p = base_params;            
    p.a2 = a2_values(i);

    fprintf('  a2 = %.3f ... ', p.a2);

    [ok, info] = assumption31_holds(p);
    if ~ok
        fprintf('infeasible, check');
        continue;
    end

    check_params_valid(p);
    [y_star, v_qvi_eval, ratio] = compute_one_sensitivity(p, settings, y_eval);

    results_a2.y_star(i) = y_star;
    results_a2.v_qvi_eval(i) = v_qvi_eval;

    results_a2.Delta_v_ratio(i) = ratio;

    fprintf('done: y*=%.4f, ratio=%.4f\n', y_star, ratio);
end

%% sensitivity to a3
% a3_values = [0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34];
a3_values = [0.24, 0.27, 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.48];

results_a3 = table();
results_a3.a3 = a3_values(:);
results_a3.y_star = NaN(numel(a3_values),1);
results_a3.Delta_v_ratio = NaN(numel(a3_values),1);

for i = 1:numel(a3_values)
    p = base_params;            
    p.a3 = a3_values(i);

    fprintf('  a3 = %.3f ... ', p.a3);

    [ok, info] = assumption31_holds(p);
    if ~ok
        fprintf('infeasible, check');
        continue;
    end

    check_params_valid(p);
    [y_star, v_qvi_eval, ratio] = compute_one_sensitivity(p, settings, y_eval);

    results_a3.y_star(i) = y_star;
    results_a3.v_qvi_eval(i) = v_qvi_eval, ;
    results_a3.Delta_v_ratio(i) = ratio;

    fprintf('done: y*=%.4f, ratio=%.4f\n', y_star, ratio);
end

%% output
fprintf('\nSensitivity to a1:\n');
disp(results_a1);

fprintf('\nSensitivity to a2:\n');
disp(results_a2);

fprintf('\nSensitivity to a3:\n');
disp(results_a3);

%% local functions
function [y_star, v_qvi_eval, Delta_v_ratio] = compute_one_sensitivity(p, settings, y_eval)

    % Suppress optional root-diagnostic messages during batch runs.
    evalc('sol_qvi = solve_bank_qvi(p, settings, ''QVI'');');
    evalc('sol_hjb = solve_bank_qvi(p, settings, ''HJB'');');

    y_star = sol_qvi.y_star;

    v_qvi_eval = interp1(sol_qvi.y, sol_qvi.v, y_eval, 'linear', 'extrap');
    v_hjb_eval = interp1(sol_hjb.y, sol_hjb.v, y_eval, 'linear', 'extrap');

    Delta_v = v_qvi_eval - v_hjb_eval;
    Delta_v_ratio = Delta_v / v_hjb_eval;
end

function [tf, info] = assumption31_holds(p)
    % Assumption 3.1 finite-value gate:
    % rho > max(mu_L, r + (mu-r)^+ * pi_inf),
    % where pi_inf is the asymptotic investment cap.
    %
    % In the unrestricted case:
    % pi_inf = 1 / max(a1,a3).
    %
    % If pi_cap is imposed:
    % pi_inf = min(1/max(a1,a3), pi_cap).

    a_bar = max(p.a1, p.a3);

    if isfield(p, 'pi_cap') && ~isempty(p.pi_cap)
        pi_inf = min(1/a_bar, p.pi_cap);
    else
        pi_inf = 1/a_bar;
    end

    growth_cap = p.r + max(p.mu - p.r, 0) * pi_inf;
    threshold = max(p.mu_L, growth_cap);
    margin = p.rho - threshold;

    tf = margin > 1e-10;

    info = struct();
    info.a_bar = a_bar;
    info.pi_inf = pi_inf;
    info.growth_cap = growth_cap;
    info.threshold = threshold;
    info.margin = margin;
end