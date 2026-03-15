clear all; clc; close all;
%% parameters
params.r       = 0.01;      
params.mu      = 0.04;      
params.mu_L    = 0.03;      
params.rho     = 0.12;      
params.sigma   = 0.08;      
params.sigma_L = 0.03;      
params.c       = 0.20;      
params.gamma   = 0.01;      
% Basel parameters
params.a1 = 0.045; %solvency
params.a2 = 0.05; %LCR
params.a3 = 0.30; %HQLA
% issuance cost
params.kappa   = 0.01;      
params.kappa_p = 0.02;      
params.pi_cap  = [];        % Empty = no cap, we can put cap = 1
check_params_valid(params);
%% discretization
settings.y_max        = 2.5;
settings.N            = 301;
settings.dt           = 0.05;
settings.max_iter     = 1500;
settings.tol          = 1e-6;
settings.pi_grid_size = 151;
%% main calculation
%compute bounds
bounds = compute_value_bounds(params);
% solve VI and 
sol_qvi = solve_bank_qvi(params, settings, 'QVI');
sol_hjb = solve_bank_qvi(params, settings, 'HJB');
% compute metrics
metrics = compute_metrics(sol_qvi, sol_hjb);
%print_metrics(metrics);
%% plots
y = sol_qvi.y;
upper_bound = y + bounds.C / bounds.rho_L;
lower_bound = y - 1;
figure(1)
% Plot 1:value Function v(y)
plot(y, sol_qvi.v, 'b-', 'LineWidth', 2, 'DisplayName', 'v_{VI}(y)');
hold on;
plot(y, sol_hjb.v, 'r--', 'LineWidth', 2, 'DisplayName', 'v_{HJB}(y)');
plot(y, lower_bound, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Lower bound: y-1');
plot(y, upper_bound, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Upper bound');
xline(sol_qvi.y_star, 'b-', 'LineWidth', 1, 'DisplayName', sprintf('y^*=%.3f', sol_qvi.y_star));
%xline(sol_hjb.y_star, 'r-', 'LineWidth', 1, 'DisplayName', sprintf('y^*_{HJB}=%.3f', sol_hjb.y_star));
xlabel('y = X/L');
ylabel('v(y)');
title('Value Function v(y)');
legend('Location', 'northwest');
xlim([1, settings.y_max]);
% remember to set reasonable y limits based on VI solution
ylim([0, max(sol_qvi.v) * 1.2]); 
grid on;
%
% Plot 2: \Delta v = v_VI - v_HJB, issuance value
figure(2)
Delta_v = sol_qvi.v - sol_hjb.v;
plot(sol_qvi.y, Delta_v, 'b-', 'LineWidth', 2);
hold on;
xline(sol_qvi.y_star, 'b--', 'LineWidth', 1); hold off;
xlabel('y = X/L');
ylabel('\Delta v = v_{VI} - v_{HJB}');
title('Value of issuance \Delta v(y)');
grid on;
%
% Plot 3: v(X,L) 
figure(3)
L_grid = linspace(0.5, 2, 50);
X_grid = linspace(0.5, 3, 50);
[X_mesh, L_mesh] = meshgrid(X_grid, L_grid);
Y_mesh = X_mesh ./ L_mesh;
V_XL = zeros(size(Y_mesh));
y_star_val = sol_qvi.y_star;
v_at_ystar = interp1(sol_qvi.y, sol_qvi.v, y_star_val, 'linear', 'extrap');
L_line = linspace(0.5, 2, 30);
v_at_1 = sol_qvi.v(1);

for i = 1:size(Y_mesh, 1)
    for j = 1:size(Y_mesh, 2)
        y_val = Y_mesh(i,j);
        % 1 <= y <= y* show only in the continuation region
        if y_val >= 1 && y_val <= y_star_val
            v_y = interp1(sol_qvi.y, sol_qvi.v, y_val, 'linear', 'extrap');
            V_XL(i,j) = L_mesh(i,j) * v_y; 
        else
            V_XL(i,j) = NaN;
        end
    end
end
contourf(X_mesh, L_mesh, V_XL, 20);
hold on;
plot(L_line, L_line, 'r-', 'LineWidth', 2);  % y=1 line
plot(y_star_val * L_line, L_line, 'g-', 'LineWidth', 2);  % y=y* line
xlabel('X ');
ylabel('L ');
title('v(X,L) ');
legend('', 'y=1 ', sprintf('y=y^*=%.2f', y_star_val), ...
       'Location', 'northwest');
colorbar;
%
