%% simulate ONE trajectory for illustration purpose
clear; clc; close all;
%%
% load the vi solution from a1
load('vi_solution.mat', 'results');
params = results.params;
sol_qvi = results.sol_qvi;

y      = sol_qvi.y;
v      = sol_qvi.v;
pi_star = sol_qvi.pi_star;
y_star = sol_qvi.y_star;

if isfield(sol_qvi, 'v_prime') && numel(sol_qvi.v_prime) == numel(y)
    v_prime = sol_qvi.v_prime;
else
    dy = y(2) - y(1);
    v_prime = zeros(size(v));
    v_prime(2:end-1) = (v(3:end) - v(1:end-2)) / (2*dy);
    v_prime(1)       = (v(2) - v(1)) / dy;
    v_prime(end)     = (v(end) - v(end-1)) / dy;
end
% compute y_post
target_slope = 1 / (1 - params.kappa_p);
[y_post, y_post_method] = pick_y_post(y, v, v_prime, y_star, params.kappa_p);
dy = y(2) - y(1);
y_trigger = params.underline_y;

% Is recapitalization actually optimal at the barrier? K subset {1} may be
% empty, in which case the bank LIQUIDATES at underline_y and never jumps to
% y_post. Gate the issuance jump on this instead of always recapitalizing.
Hv_bd = compute_H_operator(y, v, params);
active_issue = Hv_bd(1) > (params.underline_y - 1) + 1e-8;

%% simulation parameters
T_sim   = 50;% in years
dt     = 0.001;
n_steps = ceil(T_sim / dt);

% Initial conditions
Y0 = 1.20; 
L0 = 1.0;

sqrt_dt = sqrt(dt);
t_series = (0:n_steps)' * dt;
Y_series = zeros(n_steps+1, 1);
L_series = zeros(n_steps+1, 1);
X_series = zeros(n_steps+1, 1);

Y_series(1) = Y0;
L_series(1) = L0;
X_series(1) = Y0 * L0;

% Events
issuance_times    = [];
issuance_xi       = [];   
issuance_y_before = [];
issuance_y_after  = [];

dividend_times = [];
dividend_dZ    = [];     

n_issuances = 0;
total_xi    = 0;
total_dZ    = 0;

%% Trajectory simulation
Y = Y0;
for step = 1:n_steps
    t = step * dt;
    pi_opt = interp1(y, pi_star, Y, 'linear', 'extrap');
    pi_opt = max(0, min(compute_pi_bar(Y, params), pi_opt));
    [b, sigma2] = drift_and_sigma2(Y, pi_opt, params);
    sigmaY = sqrt(max(0, sigma2));
    dW = sqrt_dt * randn();
    Y_new = Y + b * dt + sigmaY * dW;
    % dividend reflection
    if Y_new > y_star
        dZ = Y_new - y_star;
        total_dZ = total_dZ + dZ;
        dividend_times = [dividend_times; t];
        dividend_dZ    = [dividend_dZ; dZ];
        Y_new = y_star;
    end
    %issuance jumps (only if recapitalization is optimal; otherwise liquidate)
    if Y_new <= y_trigger
        if ~active_issue
            % Liquidation at the distress boundary: terminate the trajectory.
            Y = params.underline_y;
            Y_series(step+1) = Y;
            L_series(step+1) = L0 * exp(params.mu_L * t);
            X_series(step+1) = Y * L_series(step+1);
            Y_series(step+2:end) = NaN;   % path ends here
            L_series(step+2:end) = NaN;
            X_series(step+2:end) = NaN;
            break;
        end
        n_issuances = n_issuances + 1;
        y_before = params.underline_y;
        xi = (y_post - ((1-params.kappa)*y_before + params.kappa)) / (1 - params.kappa_p);
        total_xi = total_xi + xi;
        issuance_times    = [issuance_times; t];
        issuance_xi       = [issuance_xi; xi];
        issuance_y_before = [issuance_y_before; y_before];
        issuance_y_after  = [issuance_y_after; y_post];
        Y_new = y_post;% jump to post-issuance target
        if Y_new > y_star
            dZ = Y_new - y_star;
            total_dZ = total_dZ + dZ;
            dividend_times = [dividend_times; t];
            dividend_dZ    = [dividend_dZ; dZ];
            Y_new = y_star;
        end
    end

    Y = max(params.underline_y, Y_new);
    Y_series(step+1) = Y;

    %for illustration if needde
    L_series(step+1) = L0 * exp(params.mu_L * t);
    X_series(step+1) = Y * L_series(step+1);
end

%% Plots
% Plot 1: Y_t
figure(1)
plot(t_series, Y_series, 'b-', 'LineWidth', 1.2);
hold on;
yline(y_star, 'g--', 'LineWidth', 2);
yline(params.underline_y, 'r--', 'LineWidth', 2);
yline(y_post,'m-.', 'LineWidth', 1.5);

% marker on the issuance
for i = 1:numel(issuance_times)
    t_iss = issuance_times(i);
    y_bef = issuance_y_before(i);
    y_aft = issuance_y_after(i);
    plot([t_iss, t_iss], [y_bef, y_aft], 'r-', 'LineWidth', 2);
    plot(t_iss, y_bef, 'rv', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    plot(t_iss, y_aft, 'r^', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
end

% marker on the reflection
for i = 1:min(100, numel(dividend_times))
    plot(dividend_times(i), y_star, 'gv', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
end

xlabel('Time (years)');
ylabel('Y(t) = X(t)/L(t)');
title(sprintf('Y(t) with controls: %d issuances, total dZ=%.4g', n_issuances, total_dZ));
legend('Y(t)', sprintf('y^* = %.3f', y_star), sprintf('y_{post}=%.3f', y_post), 'Location', 'best');

xlim([0, T_sim]);
ymin = min(params.underline_y - 0.005, min(Y_series));
ymax = max([y_star + 0.15, y_post + 0.05, max(Y_series)+0.02]);
ylim([ymin, ymax]);
grid on;
% Plot 2: cumulative issuance
figure(2)
subplot(1, 2, 1);
cum_xi = zeros(n_steps + 1, 1);
for i = 1:numel(issuance_times)
    idx = find(t_series >= issuance_times(i), 1);
    if ~isempty(idx)
        cum_xi(idx:end) = cum_xi(idx:end) + issuance_xi(i);
    end
end
stairs(t_series, cum_xi, 'r-', 'LineWidth', 1.8);
xlabel('Time (years)');
ylabel('Cumulative \xi ');
title(sprintf('Total capital issued: %.4g', total_xi));
xlim([0, T_sim]);
grid on;
% Plot 3: cumulative dividends
subplot(1, 2, 2);
cum_dZ = zeros(n_steps + 1, 1);
for i = 1:numel(dividend_times)
    idx = find(t_series >= dividend_times(i), 1);
    if ~isempty(idx)
        cum_dZ(idx:end) = cum_dZ(idx:end) + dividend_dZ(i);
    end
end
plot(t_series, cum_dZ, 'g-', 'LineWidth', 1.8);
xlabel('Time (years)');
ylabel('Cumulative dZ');
title(sprintf('Total dividends: %.4g', total_dZ));
xlim([0, T_sim]);
grid on;
