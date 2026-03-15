%% simulate ONE trajectory for illustration purpose
clear; clc; close all;
%%
% load the vi solution
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
eps_hit = max(1e-6, 0.5*dy);
y_trigger = 1 + eps_hit;
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
    %issuance jumps
    if Y_new <= y_trigger
        n_issuances = n_issuances + 1;
        y_before = 1.0;
        xi = (y_post - y_before) / (1 - params.kappa_p);
        total_xi = total_xi + xi;
        issuance_times    = [issuance_times; t];
        issuance_xi       = [issuance_xi; xi];
        issuance_y_before = [issuance_y_before; y_before];
        issuance_y_after  = [issuance_y_after; y_post];
        Y_new = y_post;%always jumps to y_post
        if Y_new > y_star
            dZ = Y_new - y_star;
            total_dZ = total_dZ + dZ;
            dividend_times = [dividend_times; t];
            dividend_dZ    = [dividend_dZ; dZ];
            Y_new = y_star;
        end
    end

    Y = max(1.0, Y_new);
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
yline(1.0,   'r--', 'LineWidth', 2);
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
ymin = min(0.995, min(Y_series));
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

%% local function
function [y_post, method] = pick_y_post(y, v, v_prime, y_star, kappa_p)
    target = 1 / (1 - kappa_p);
    %y_post need to be in the continuation region
    I = find(y <= y_star);
    if numel(I) < 5
        [y_post, method] = pick_y_post_from_H_discrete(y, v, kappa_p);
        return;
    end

    yI = y(I);
    g  = v_prime(I) - target;
    i0 = min(5, numel(yI)-1);
    y_post = NaN;
    for j = i0+1:numel(yI)
        if g(j-1) > 0 && g(j) <= 0
            denom = g(j-1) - g(j);
            if abs(denom) > 1e-14
                alpha = g(j-1) / denom; 
                y_post = yI(j-1) + alpha * (yI(j) - yI(j-1));
            else
                y_post = 0.5 * (yI(j-1) + yI(j));
            end
            break;
        end
    end
    if isnan(y_post)
        warning('y_post not found');
        y_post = y_star;
        method = 'fallback';
    end
end
