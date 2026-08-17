%% issuance cost test on kappa'
clear; clc; close all;
%% baseline
base_params = struct();
base_params.r       = 0.02;
base_params.mu      = 0.04;
base_params.mu_L    = 0.03;
base_params.rho     = 0.12;
base_params.sigma   = 0.08;
base_params.sigma_L = 0.03;
base_params.c       = 0.20;
base_params.gamma   = 0.02;
base_params.kappa   = 0.01;
base_params.regime  = 'both';
base_params.pi_cap  = 1.0;       % tight ceiling
base_params.a2 = 0.05;
base_params.a3 = 0.15;

%% sweep
kappa_p_grid = [0.02, 0.05, 0.10, 0.20];
a1_grid      = [0.045, 0.060, 0.080, 0.100, 0.120, 0.150];

h_master     = 0.005;
uly_step     = 0.010;
uly_grid     = (round(1.010/uly_step) : round(1.180/uly_step)) * uly_step;

y0        = 1.20;
T_horizon = 5.0;

occ_levels = [1.05, 1.10, 1.15];

y_max_safety = 1.15;
y_max_pad    = 0.50;

settings = struct('dt',0.05,'max_iter',3000,'tol',1e-9, ...
    'pi_grid_size',151,'verbose',false,'compute_analytical_ystar',false);

n_paths = 5000;
dt_mc   = 0.05;
seed    = 12345;

for f = {'solve_bank_qvi','compute_pi_bar','drift_and_sigma2','assumption31_holds'}
    if isempty(which(f{1})), error('%s not on the path.', f{1}); end
end

%% enumerate
W = {};
for ik = 1:numel(kappa_p_grid)
    for i1 = 1:numel(a1_grid)
        for iu = 1:numel(uly_grid)
            p = base_params;
            p.kappa_p     = kappa_p_grid(ik);
            p.a1          = a1_grid(i1);
            p.underline_y = uly_grid(iu);

            [ok, info] = assumption31_holds(p);
            if ~ok, continue; end

            y_target = max(y_max_safety*info.y_star_bound, p.underline_y + y_max_pad);
            y_max = p.underline_y + ceil((y_target - p.underline_y)/h_master)*h_master;

            s = settings;
            s.y_max = y_max;
            s.N     = round((y_max - p.underline_y)/h_master) + 1;

            W{end+1} = struct('p',p,'s',s); %#ok<SAGROW>
        end
    end
end

nW = numel(W);
fprintf('Solves: %d  (N = %d)\n\n', nW, W{1}.s.N);

%% solve
res = cell(nW,1); okv = false(nW,1);

t0 = tic;
for w = 1:nW
    if okv(w), continue; end
    try
        res{w} = solve_and_measure(W{w}.p, W{w}.s, y0, T_horizon, ...
            n_paths, dt_mc, seed, occ_levels);
        okv(w) = true;
    catch ME
        fprintf('  %d FAILED: %s\n', w, ME.message);
    end
    if mod(w, 40) == 0
        el = toc(t0);
        fprintf('  %4d/%4d | %5.1f min | ~%.0f min left\n', ...
            w, nW, el/60, el/w*(nW-w)/60);
        save(checkpoint, 'res', 'okv', '-v7.3');
    end
end
save(checkpoint, 'res', 'okv', '-v7.3');

rows = {};
for w = 1:nW
    if ~okv(w), continue; end
    r = res{w}; p = W{w}.p;
    rows{end+1} = [p.kappa_p, p.a1, p.underline_y, r.v_y0, r.y_star, ...
        r.y_star - p.underline_y, r.y_post, r.xi_post, double(r.active_issue), ...
        r.E_N, r.E_issue, r.mean_Y, r.occ(:)']; %#ok<SAGROW>
end
M = cell2mat(rows(:));

vn = [{'kappa_p','a1','underline_y','v_y0','y_star','band','y_post','xi_post', ...
       'active_issue','E_N','E_disc_issuance','mean_Y'}, ...
      arrayfun(@(L) sprintf('occ_below_%03d', round(1000*L)), occ_levels, ...
               'UniformOutput', false)];
T = array2table(M, 'VariableNames', vn);

%% report
fprintf('\n%s\nbanks response with recapitalization costs\n%s\n', ...
    repmat('=',1,78), repmat('=',1,78));
fprintf('  kappa''   a1      y^B    band    E[N]    E[iss]   occ<%.2f   v(y0)\n', ...
    occ_levels(2));
occ_col = sprintf('occ_below_%03d', round(1000*occ_levels(2)));
for ik = 1:numel(kappa_p_grid)
    for i1 = 1:numel(a1_grid)
        S = T(abs(T.kappa_p-kappa_p_grid(ik))<1e-12 & abs(T.a1-a1_grid(i1))<1e-12, :);
        if isempty(S), continue; end
        [~, j] = max(S.v_y0);
        fprintf('  %.3f  %.3f  %.3f  %.4f  %6.2f  %7.4f  %8.4f  %.5f\n', ...
            kappa_p_grid(ik), a1_grid(i1), S.underline_y(j), S.band(j), ...
            S.E_N(j), S.E_disc_issuance(j), S.(occ_col)(j), S.v_y0(j));
    end
    fprintf('\n');
end
fprintf('\nSaved %s\n', csv_file);

%% local functions

function res = solve_and_measure(p, s, y0, T_horizon, n_paths, dt, seed, occ_levels)

    evalc('sol = solve_bank_qvi(p, s, ''QVI'');');

    [d0, j0] = min(abs(sol.y - y0));
    if d0 > 1e-9, error('y0 not a grid node.'); end

    res = struct();
    res.v_y0   = sol.v(j0);
    res.y_star = sol.y_star;

    [y_post, xi_post, active] = boundary_post_target(sol.y, sol.v, p);
    res.y_post = y_post; res.xi_post = xi_post; res.active_issue = active;

    rng(seed, 'twister');
    yg = sol.y(:); pig = sol.pi_star(:); ystar = sol.y_star;
    rho_L = p.rho - p.mu_L;

    n_steps = ceil(T_horizon/dt); sdt = sqrt(dt);
    Y = y0*ones(n_paths,1);
    nint = zeros(n_paths,1); iss = zeros(n_paths,1);
    occ = zeros(n_paths, numel(occ_levels));
    intY = zeros(n_paths,1);
    alive = true(n_paths,1);

    for step = 1:n_steps
        ia = find(alive);
        if isempty(ia), break; end
        t = step*dt;
        Ya = Y(ia);

        for iL = 1:numel(occ_levels)
            occ(ia,iL) = occ(ia,iL) + dt*(Ya < occ_levels(iL));
        end
        intY(ia) = intY(ia) + dt*Ya;

        pi_o = interp1(yg, pig, Ya, 'linear', 'extrap');
        try
            cap = compute_pi_bar(Ya, p);
            if isscalar(cap), cap = repmat(cap, size(Ya)); end
            if numel(cap) ~= numel(Ya), error('nv'); end
        catch
            cap = arrayfun(@(yy) compute_pi_bar(yy,p), Ya);
        end
        pi_o = max(0, min(cap(:), pi_o));

        try
            [b, s2] = drift_and_sigma2(Ya, pi_o, p);
            if isscalar(b),  b  = repmat(b, size(Ya)); end
            if isscalar(s2), s2 = repmat(s2, size(Ya)); end
            if numel(b) ~= numel(Ya), error('nv'); end
        catch
            b = zeros(numel(Ya),1); s2 = zeros(numel(Ya),1);
            for q = 1:numel(Ya)
                [b(q), s2(q)] = drift_and_sigma2(Ya(q), pi_o(q), p);
            end
        end
        b = b(:); s2 = max(0, s2(:));

        Yn = min(Ya + b*dt + sqrt(s2).*sdt.*randn(size(Ya)), ystar);

        U  = rand(numel(Ya),1);
        cr = (Yn <= p.underline_y);
        cand = ~cr & (Ya > p.underline_y) & (Yn > p.underline_y) & (s2 > 1e-14);
        if any(cand)
            j = find(cand);
            ex = -2*(Ya(j)-p.underline_y).*(Yn(j)-p.underline_y)./(s2(j)*dt);
            cr(j) = U(j) < exp(min(0, ex));
        end

        if any(cr)
            jc = ia(cr);
            nint(jc) = nint(jc) + 1;
            if active
                iss(jc) = iss(jc) + exp(-rho_L*t)*xi_post;
                Yn(cr)  = y_post;
            else
                alive(jc) = false;
                Yn(cr) = p.underline_y;
            end
        end
        Y(ia) = Yn;
    end

    res.E_N     = mean(nint);
    res.E_issue = mean(iss);
    res.occ     = mean(occ,1) / T_horizon;   % fraction of the horizon
    res.mean_Y  = mean(intY) / T_horizon;
end

function [y_post, xi_post, active] = boundary_post_target(y, v, p)
    y_lo = p.underline_y; om = 1 - p.kappa_p;
    obj = v - (y - y_lo)/om - (p.kappa*p.kappa_p/om)*(y_lo - 1);
    obj(y <= y_lo) = -Inf;
    [Hval, j] = max(obj);
    active = Hval > (y_lo - 1) + 1e-8;
    if active
        y_post  = y(j);
        xi_post = max((y_post - (1-p.kappa)*y_lo - p.kappa)/om, 0);
    else
        y_post = NaN; xi_post = NaN;
    end
end
