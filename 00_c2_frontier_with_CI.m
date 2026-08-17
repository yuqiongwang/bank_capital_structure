%% Pareto frontier with Monte Carlo error bars and CI
% the version without error bar is actually present in the paper now
% but we do this one for the table to disclose the CI
clear; clc; close all;
%% parameters (baseline, r > r_L)
params.r = 0.02; params.mu = 0.04; params.mu_L = 0.03;
params.sigma = 0.08; params.sigma_L = 0.03; params.c = 0.20;
params.gamma = 0.02; params.rho = 0.12;
params.kappa = 0.01; params.kappa_p = 0.02;
params.pi_cap = [];          
params.underline_y = 1.02;
params.a1 = 0.045; params.a2 = 0.05; params.a3 = 0.25;

a1_grid = linspace(0.045,0.12,20);
% don't use the original a1 grid, let us refine the grid 
%a1_grid = [0.045, 0.06, 0.08, 0.10, 0.12];
a2_grid = [0.05, 0.10];
a3_grid = [0.15, 0.25, 0.50];

T_horizon = 5.0;  y0 = 1.20;
settings = struct('N',301,'y_max',2.5,'dt',0.05,'max_iter',1500,'tol',1e-6,'pi_grid_size',151);
check_params_valid(params);

n_paths   = 2000;  
dt        = 0.05;
z_score   = 1.96;    % 95% CI
mc_seed   = 20240517;   

n1=numel(a1_grid); n2=numel(a2_grid); n3=numel(a3_grid); n_total=n1*n2*n3;
R = struct();
[R.a1,R.a2,R.a3,R.y_star,R.v_y0,R.P,R.P_lo,R.P_hi] = deal(zeros(n_total,1));

idx = 0;
for i1=1:n1
 for i2=1:n2
  for i3=1:n3
    idx = idx+1;
    params = params;
    params.a1=a1_grid(i1); params.a2=a2_grid(i2); params.a3=a3_grid(i3);
    [ok_feas, ~] = check_params_valid(params);
    if ~ok_feas
        R.a1(idx)=params.a1; R.a2(idx)=params.a2; R.a3(idx)=params.a3;
        [R.y_star(idx),R.v_y0(idx),R.P(idx),R.P_lo(idx),R.P_hi(idx)]=deal(NaN);
        continue;
    end
    try
        sol  = solve_bank_qvi(params, settings, 'QVI');
        v_y0 = interp1(sol.y, sol.v, y0, 'linear', 'extrap');
        rng(mc_seed);   
        Phat = estimate_survival_prob(sol, params, y0, T_horizon, n_paths, dt);
        x  = round(Phat * n_paths);
        [plo, phi] = wilson_ci(x, n_paths, z_score);
        R.a1(idx)=params.a1; R.a2(idx)=params.a2; R.a3(idx)=params.a3;
        R.y_star(idx)=sol.y_star; R.v_y0(idx)=v_y0;
        R.P(idx)=Phat; R.P_lo(idx)=plo; R.P_hi(idx)=phi;
    catch
        R.a1(idx)=params.a1; R.a2(idx)=params.a2; R.a3(idx)=params.a3;
        [R.y_star(idx),R.v_y0(idx),R.P(idx),R.P_lo(idx),R.P_hi(idx)]=deal(NaN);
    end
  end
 end
end

valid = ~isnan(R.v_y0) & ~isnan(R.P);

%% Point-estimate Pareto frontier 
is_par = pareto_mask(R.v_y0, R.P, valid);

%% Robust frontier: a point is "robustly efficient" if no other point
% dominates it even when we credit its upper survival bound
% and penalize with its lowe survival bound r_l
is_robust = robust_pareto_mask(R.v_y0, R.P_lo, R.P_hi, valid);

%% report
T = table(R.a1, R.a2, R.a3, R.y_star, R.v_y0, R.P, R.P_lo, R.P_hi, is_par, is_robust, ...
    'VariableNames', {'a1','a2','a3','y_star','v_y0','P','P_lo','P_hi','pareto','robust'});
T = sortrows(T(valid,:), 'P');
fprintf('\nFrontier with %d MC paths, 95%% CIs:\n', n_paths);
disp(T);
fprintf('Point-estimate frontier: %d points;  robust frontier: %d points.\n', ...
    sum(is_par), sum(is_robust));

%% plot with horizontal error bars
figure(1); hold on;
errorbar(R.P(valid), R.v_y0(valid), [], [], ...
    R.P(valid)-R.P_lo(valid), R.P_hi(valid)-R.P(valid), ...
    'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4, 'CapSize', 0);
scatter(R.P(is_par),    R.v_y0(is_par),    70, 'b', 'filled', 'DisplayName','point frontier');
scatter(R.P(is_robust), R.v_y0(is_robust), 30, 'r', 'filled', 'DisplayName','robust frontier');
xlabel('P(\tau \geq T)  (95% CI)'); ylabel('v(y_0)');
title('Profit-safety frontier with CI');
legend('all configs with 95\% CI', 'point frontier','robust frontier','Location','southwest'); grid on;

%% ---- local functions ----
function [lo,hi] = wilson_ci(x, n, z)
    phat = x/n;
    denom = 1 + z^2/n;
    centre = (phat + z^2/(2*n)) / denom;
    halfw  = (z*sqrt(phat*(1-phat)/n + z^2/(4*n^2))) / denom;
    lo = max(0, centre - halfw);
    hi = min(1, centre + halfw);
end

function mask = pareto_mask(v, P, valid)
    n = numel(v); mask = false(n,1);
    vi = find(valid);
    for ii = 1:numel(vi)
        i = vi(ii); dom = false;
        for jj = 1:numel(vi)
            j = vi(jj);
            if j~=i && v(j)>=v(i) && P(j)>=P(i) && (v(j)>v(i) || P(j)>P(i))
                dom = true; break;
            end
        end
        mask(i) = ~dom;
    end
end

function mask = robust_pareto_mask(v, P_lo, P_hi, valid)
    n = numel(v); mask = false(n,1);
    vi = find(valid);
    for ii = 1:numel(vi)
        i = vi(ii); dom = false;
        for jj = 1:numel(vi)
            j = vi(jj);
            if j~=i && v(j)>=v(i) && P_lo(j) >= P_hi(i) && ...
               (v(j)>v(i) || P_lo(j)>P_hi(i))
                dom = true; break;
            end
        end
        mask(i) = ~dom;
    end
end
