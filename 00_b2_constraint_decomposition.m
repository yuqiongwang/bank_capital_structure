%% constraint decomposition
% Compare four regulatory regimes: none, lcr, solvency, both
% the results must be nested. double check
clear; clc; close all;
%% parameters
base_params = struct();

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
base_params.underline_y = 1.02;

pi_cap_common = [];
base_params.pi_cap = pi_cap_common;

%% discretization
settings = struct();
settings.y_max        = 2.5;
settings.N            = 301;
settings.dt           = 0.05;
settings.max_iter     = 1500;
settings.tol          = 1e-6;
settings.pi_grid_size = 151;
settings.verbose      = false;
settings.compute_analytical_ystar = false;

%% evaluation point and regimes
y_eval  = 1.20;
regimes = {'none','solvency','lcr','both'};
nR      = numel(regimes);

%% plotting styles 
% to make sure we can also distinguish when the print is not in color 
plot_colors = lines(nR);
line_styles = {'-', '--', '-.', ':'};
markers     = {'o', 's', '^', 'd'};

line_width  = 1.8;
marker_size = 6;
n_markers   = 10;

%% diagnostics: basic conditions
r_L = base_params.mu_L - base_params.gamma;
if base_params.r <= r_L
    warning('r > r_L violated: r=%.4f, r_L=%.4f.', ...
        base_params.r, r_L);
else
    fprintf('condition r > r_L: ok (r=%.4f, r_L=%.4f)\n', ...
        base_params.r, r_L);
end

check_params_valid(base_params);
fprintf('Using compute_pi_bar at: %s\n', which('compute_pi_bar'));
src = fileread(which('compute_pi_bar'));
%% check caps at one point
y_check = 1.50;
pcheck = zeros(nR,1);

for k = 1:nR
    p = base_params;
    p.regime = regimes{k};
    pcheck(k) = compute_pi_bar(y_check, p);
end

fprintf(['\npi_bar(%.2f) by regime [none solvency lcr both] = ', ...
         '[%.4f %.4f %.4f %.4f]\n'], ...
    y_check, pcheck(1), pcheck(2), pcheck(3), pcheck(4));

%% cap nesting diagnostic
ytest = linspace(base_params.underline_y, settings.y_max, 200);
caps = zeros(numel(ytest), nR);

for k = 1:nR
    p = base_params;
    p.regime = regimes{k};
    caps(:,k) = arrayfun(@(yy) compute_pi_bar(yy, p), ytest);
end

fprintf('all these should be <=0 up to :\n');
fprintf('max(sol - none) = %.4e\n', max(caps(:,2) - caps(:,1)));
fprintf('max(lcr - none) = %.4e\n', max(caps(:,3) - caps(:,1)));
fprintf('max(both - sol) = %.4e\n', max(caps(:,4) - caps(:,2)));
fprintf('max(both - lcr) = %.4e\n', max(caps(:,4) - caps(:,3)));

if max(caps(:,2) - caps(:,1)) > 1e-10 || ...
   max(caps(:,3) - caps(:,1)) > 1e-10 || ...
   max(caps(:,4) - caps(:,2)) > 1e-10 || ...
   max(caps(:,4) - caps(:,3)) > 1e-10
    warning('the checking does not pass, check gain the parameters.');
end

%% solve each regime
res = struct();
res.regime        = regimes(:);
res.y_star        = NaN(nR,1);
res.y_post        = NaN(nR,1);
res.xi_post       = NaN(nR,1);
res.v_at_eval     = NaN(nR,1);
res.v_hjb_eval    = NaN(nR,1);
res.G             = NaN(nR,1);

sols_qvi = cell(nR,1);
sols_hjb = cell(nR,1);

for k = 1:nR
    p = base_params;
    p.regime = regimes{k};

    fprintf('\nregime = %-9s ... ', p.regime);

    try
        % Suppress optional root-diagnostic messages.
        evalc('sol_qvi = solve_bank_qvi(p, settings, ''QVI'');');
        evalc('sol_hjb = solve_bank_qvi(p, settings, ''HJB'');');

        sols_qvi{k} = sol_qvi;
        sols_hjb{k} = sol_hjb;

        v_qvi = interp1(sol_qvi.y, sol_qvi.v, ...
            y_eval, 'linear', 'extrap');
        v_hjb = interp1(sol_hjb.y, sol_hjb.v, ...
            y_eval, 'linear', 'extrap');

        [y_post, xi_post, active_issue] = ...
            get_boundary_post_target(sol_qvi.y, sol_qvi.v, p);

        res.y_star(k)     = sol_qvi.y_star;
        res.y_post(k)     = y_post;
        res.xi_post(k)    = xi_post;
        res.v_at_eval(k)  = v_qvi;
        res.v_hjb_eval(k) = v_hjb;
        res.G(k)          = (v_qvi - v_hjb) / v_hjb;

        fprintf(['done: y*=%.4f, y_post=%.4f, xi=%.4f, ', ...
                 'v(%.2f)=%.4f, G=%.4f\n'], ...
            sol_qvi.y_star, y_post, xi_post, ...
            y_eval, v_qvi, res.G(k));

    catch ME
        fprintf('failed: %s\n', ME.message);
    end
end

%% table
T = table(res.regime, res.y_star, res.y_post, res.xi_post, ...
          res.v_at_eval, res.v_hjb_eval, res.G, ...
    'VariableNames', {'regime','y_star','y_post','xi_post', ...
                      'v_QVI_y0','v_HJB_y0','G_y0'});

fprintf('\nConstraint decomposition at y0 = %.2f:\n', y_eval);
disp(T);

%% value nesting diagnostic
fprintf('\nValue nesting diagnostics:\n');

idx_none = find(strcmp(regimes, 'none'), 1);
idx_sol  = find(strcmp(regimes, 'solvency'), 1);
idx_lcr  = find(strcmp(regimes, 'lcr'), 1);
idx_both = find(strcmp(regimes, 'both'), 1);

if ~isempty(sols_qvi{idx_none})
    v_none = sols_qvi{idx_none}.v;

    for k = 2:nR
        if ~isempty(sols_qvi{k})
            diff_val = v_none - sols_qvi{k}.v;
            fprintf('min(v_none - v_%s) = %.4e\n', ...
                regimes{k}, min(diff_val));
        end
    end

    if ~isempty(sols_qvi{idx_both}) && ...
       ~isempty(sols_qvi{idx_sol})
        fprintf('min(v_solvency - v_both) = %.4e\n', ...
            min(sols_qvi{idx_sol}.v - sols_qvi{idx_both}.v));
    end

    if ~isempty(sols_qvi{idx_both}) && ...
       ~isempty(sols_qvi{idx_lcr})
        fprintf('min(v_lcr - v_both) = %.4e\n', ...
            min(sols_qvi{idx_lcr}.v - sols_qvi{idx_both}.v));
    end
end


% %% plot: value functions by regime
% figure(1);
% clf;
% set(gcf, 'Color', 'w');
% hold on;
% 
% for k = 1:nR
%     if ~isempty(sols_qvi{k})
%         plot_series_with_markers( ...
%             sols_qvi{k}.y, sols_qvi{k}.v, ...
%             plot_colors(k,:), ...
%             line_styles{k}, ...
%             markers{k}, ...
%             line_width, ...
%             marker_size, ...
%             n_markers, ...
%             regimes{k});
%     end
% end
% 
% xlabel('$y=X/L$', 'Interpreter', 'latex');
% ylabel('$v(y)$', 'Interpreter', 'latex');
% title('Value function by regulatory regime');
% legend('Location', 'northwest', 'Interpreter', 'none');
% grid on;
% box on;
% xlim([base_params.underline_y, settings.y_max]);

%% plot: zoomed value functions
figure(1);
clf;
set(gcf, 'Color', 'w');
hold on;

for k = 1:nR
    if ~isempty(sols_qvi{k})
        plot_series_with_markers( ...
            sols_qvi{k}.y, sols_qvi{k}.v, ...
            plot_colors(k,:), ...
            line_styles{k}, ...
            markers{k}, ...
            line_width, ...
            marker_size, ...
            n_markers, ...
            regimes{k});
    end
end

xlabel('$y=X/L$', 'Interpreter', 'latex');
ylabel('$v(y)$', 'Interpreter', 'latex');
title('Value function by regulatory regime');% before adding to paper
% manually change the scheme in the title
legend('Location', 'northwest', 'Interpreter', 'none');
grid on;
box on;
xlim([base_params.underline_y + 0.02, 1.80]);
%% local functions
function [y_post, xi_post, active_issue] = ...
        get_boundary_post_target(y, v, p)

    y_lo = p.underline_y;
    om = 1 - p.kappa_p;

    obj = v ...
        - (y - y_lo) / om ...
        - (p.kappa * p.kappa_p / om) * (y_lo - 1);

    obj(y <= y_lo) = -Inf;

    [Hval, j] = max(obj);

    liq = y_lo - 1;
    active_issue = Hval > liq + 1e-8;

    if active_issue
        y_post = y(j);
        xi_post = ...
            (y_post - (1 - p.kappa) * y_lo - p.kappa) / om;
        xi_post = max(xi_post, 0);
    else
        y_post = NaN;
        xi_post = NaN;
    end
end

function plot_series_with_markers( ...
        x, y, color_value, line_style, marker_style, ...
        line_width, marker_size, n_markers, display_name)

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if isempty(x)
        return;
    end

    marker_idx = unique(round(linspace(1, numel(x), ...
        min(n_markers, numel(x)))));
    plot(x, y, ...
        'Color', color_value, ...
        'LineStyle', line_style, ...
        'LineWidth', line_width, ...
        'HandleVisibility', 'off');

    % markers.
    plot(x(marker_idx), y(marker_idx), ...
        'Color', color_value, ...
        'LineStyle', 'none', ...
        'Marker', marker_style, ...
        'MarkerSize', marker_size, ...
        'MarkerFaceColor', 'none', ...
        'HandleVisibility', 'off');

    % legends
        plot(NaN, NaN, ...
        'Color', color_value, ...
        'LineStyle', line_style, ...
        'LineWidth', line_width, ...
        'Marker', marker_style, ...
        'MarkerSize', marker_size, ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', display_name);
end
