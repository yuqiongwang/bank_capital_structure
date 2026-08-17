%% issuance cost test on kappa
% similar to d2 but for kappa
clear; clc; close all;
%% execution controls
use_parallel  = true;     
chunk_size    = 60;
force_restart = false;

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

base_params.kappa   = 0.01;   
base_params.kappa_p = 0.02;   % fixed

base_params.regime = 'both';
base_params.pi_cap = 1.0;    
base_params.a2 = 0.05;
base_params.a3 = 0.15;

%% grids
kappa_grid = [0.00, 0.01, 0.05, 0.10, 0.20];

a1_grid = [0.045, 0.060, 0.080, 0.100, 0.120, 0.150];

h_master = 0.005;
uly_step = 0.010;
uly_grid = (round(1.010/uly_step):round(1.180/uly_step))*uly_step;


y0        = 1.20;
T_horizon = 5.0;

y_max_safety = 1.15;
y_max_pad    = 0.50;

settings = struct('dt',0.05,'max_iter',3000,'tol',1e-9, ...
    'pi_grid_size',151,'verbose',false,'compute_analytical_ystar',false);

n_paths = 5000;
dt_mc   = 0.05;
seed    = 12345;

res = cell(nW,1);
okv = false(nW,1);

%% enumerate QVI work
W = {};
for ik = 1:numel(kappa_grid)
    for i1 = 1:numel(a1_grid)
        for iu = 1:numel(uly_grid)
            p = base_params;
            p.kappa       = kappa_grid(ik);
            p.a1          = a1_grid(i1);
            p.underline_y = uly_grid(iu);

            [ok, info] = assumption31_holds(p);
            if ~ok
                continue;
            end

            y_target = max(y_max_safety*info.y_star_bound, ...
                           p.underline_y + y_max_pad);
            y_max = p.underline_y + ...
                ceil((y_target-p.underline_y)/h_master)*h_master;

            s = settings;
            s.y_max = y_max;
            s.N     = round((y_max-p.underline_y)/h_master) + 1;

            W{end+1} = struct('p',p,'s',s); %#ok<SAGROW>
        end
    end
end

nW = numel(W);
fprintf('kappa result: %d VI solves.\n', nW);
fprintf('  kappa values: %s\n', mat2str(kappa_grid));
fprintf('  a1 values:    %s\n', mat2str(a1_grid));
fprintf('  trigger grid: %.3f to %.3f (%d points)\n\n', ...
    min(uly_grid), max(uly_grid), numel(uly_grid));

%% VI solve
has_pct = license('test','Distrib_Computing_Toolbox');
do_parallel = use_parallel && has_pct;

if do_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool;
        catch ME
            warning('cannot start', ME.message);
            do_parallel = false;
        end
    end
end

todo = find(~okv);
t0 = tic;

for ib = 1:chunk_size:numel(todo)
    idx = todo(ib:min(ib+chunk_size-1,numel(todo)));
    tmp = cell(numel(idx),1);
    tmp_ok = false(numel(idx),1);
    tmp_msg = strings(numel(idx),1);

    if do_parallel
        parfor jj = 1:numel(idx)
            w = idx(jj);
            try
                tmp{jj} = solve_qvi_only(W{w}.p, W{w}.s, y0);
                tmp_ok(jj) = true;
                tmp_msg(jj) = "ok";
            catch ME
                tmp_msg(jj) = string(ME.message);
            end
        end
    else
        for jj = 1:numel(idx)
            w = idx(jj);
            try
                tmp{jj} = solve_qvi_only(W{w}.p, W{w}.s, y0);
                tmp_ok(jj) = true;
                tmp_msg(jj) = "ok";
            catch ME
                tmp_msg(jj) = string(ME.message);
            end
        end
    end

    for jj = 1:numel(idx)
        w = idx(jj);
        if tmp_ok(jj)
            res{w} = tmp{jj};
            okv(w) = true;
        else
            fprintf('  solve %d failed: %s\n', w, tmp_msg(jj));
        end
    end

    save(checkpoint, 'res', 'okv', 'signature', '-v7.3');

    done = nnz(okv);
    elapsed = toc(t0);
    fprintf('  VI %4d/%4d done | elapsed %.1f min\n', ...
        done, nW, elapsed/60);
end

%%
rows = {};
work_index = [];
for w = 1:nW
    if ~okv(w), continue; end
    p = W{w}.p;
    r = res{w};

    rows{end+1} = [ ...
        p.kappa, p.kappa_p, p.a1, p.underline_y, ...
        r.v_y0, r.y_star, r.band, r.y_post, r.xi_post, ...
        double(r.active_issue)]; 
    work_index(end+1,1) = w;
end

M = cell2mat(rows(:));
Full = array2table(M, 'VariableNames', { ...
    'kappa','kappa_p','a1','underline_y','v_y0','y_star','band', ...
    'y_post','xi_post','active_issue'});
Full.work_index = work_index;

writetable(Full, full_csv);

%% select exact bank best responses 
BR = table();

for ik = 1:numel(kappa_grid)
    for i1 = 1:numel(a1_grid)
        S = Full(abs(Full.kappa-kappa_grid(ik))<1e-12 & ...
                 abs(Full.a1-a1_grid(i1))<1e-12, :);
        if isempty(S), continue; end

        S = sortrows(S,'underline_y');
        vmax = max(S.v_y0);

        jj = find(abs(S.v_y0-vmax) <= 1e-13*max(1,abs(vmax)));
        j = jj(1);

        r = S(j,:);
        r.profile_spread_bp = 1e4*(max(S.v_y0)-min(S.v_y0));
        r.at_grid_edge = ...
            abs(r.underline_y-min(S.underline_y))<1e-12 || ...
            abs(r.underline_y-max(S.underline_y))<1e-12;

        jlow = find(abs(S.underline_y-1.010)<1e-12,1);
        if isempty(jlow)
            r.trigger_gain_vs_1p01_bp = NaN;
        else
            r.trigger_gain_vs_1p01_bp = 1e4*(r.v_y0-S.v_y0(jlow));
        end

        BR = [BR; r]; 
    end
end

%% simulate the best response
BR.E_N = NaN(height(BR),1);
BR.E_disc_issuance = NaN(height(BR),1);

fprintf('\nSimulating %d best-response policies only...\n', height(BR));
t_mc = tic;

for j = 1:height(BR)
    w = BR.work_index(j);
    p = W{w}.p;
    r = res{w};

    try
        mc = simulate_policy(r, p, y0, T_horizon, n_paths, dt_mc, seed);
        BR.E_N(j) = mc.E_N;
        BR.E_disc_issuance(j) = mc.E_issue;
    catch ME
        fprintf('  MC row %d FAILED: %s\n', j, ME.message);
    end

    if mod(j,6)==0 || j==height(BR)
        fprintf('  MC %2d/%2d | elapsed %.1f min\n', ...
            j, height(BR), toc(t_mc)/60);
    end
end

writetable(BR, br_csv);

%% compact summary by kappa
Summary = table();
for ik = 1:numel(kappa_grid)
    S = BR(abs(BR.kappa-kappa_grid(ik))<1e-12,:);
    if isempty(S), continue; end

    row = table();
    row.kappa = kappa_grid(ik);
    row.median_underline_y_B = median(S.underline_y,'omitnan');
    row.min_underline_y_B    = min(S.underline_y);
    row.max_underline_y_B    = max(S.underline_y);
    row.median_E_N           = median(S.E_N,'omitnan');
    row.median_xi_post       = median(S.xi_post,'omitnan');
    row.median_E_disc_issuance = median(S.E_disc_issuance,'omitnan');
    row.all_at_low_trigger = all(abs(S.underline_y-1.010)<1e-12);
    row.any_at_grid_edge = any(S.at_grid_edge);

    Summary = [Summary; row]; %#ok<AGROW>
end

writetable(Summary, summary_csv);

%% local functions

function r = solve_qvi_only(p, s, y0)
    evalc('sol = solve_bank_qvi(p, s, ''QVI'');');
    [d0,j0] = min(abs(sol.y-y0));
    if d0 > 1e-9
        error('y0=%.6f is not a grid node; nearest distance %.3g.', y0, d0);
    end

    [y_post, xi_post, active] = boundary_post_target(sol.y, sol.v, p);

    r = struct();
    r.v_y0 = sol.v(j0);
    r.y_star = sol.y_star;
    r.band = sol.y_star-p.underline_y;
    r.y_post = y_post;
    r.xi_post = xi_post;
    r.active_issue = active;

    % Retain only what is required for MC at the selected policies.
    r.y_grid = sol.y(:);
    r.pi_star = sol.pi_star(:);
end

function mc = simulate_policy(r, p, y0, T_horizon, n_paths, dt, seed)

    rng(seed,'twister');

    yg = r.y_grid(:);
    pig = r.pi_star(:);
    ystar = r.y_star;
    rho_L = p.rho-p.mu_L;

    n_steps = ceil(T_horizon/dt);
    sdt = sqrt(dt);

    Y = y0*ones(n_paths,1);
    nint = zeros(n_paths,1);
    iss = zeros(n_paths,1);
    alive = true(n_paths,1);

    for step = 1:n_steps
        ia = find(alive);
        if isempty(ia), break; end

        t = step*dt;
        Ya = Y(ia);

        pi_o = interp1(yg,pig,Ya,'linear','extrap');

        try
            cap = compute_pi_bar(Ya,p);
            if isscalar(cap), cap = repmat(cap,size(Ya)); end
            if numel(cap) ~= numel(Ya), error('not vectorized'); end
        catch
            cap = arrayfun(@(yy) compute_pi_bar(yy,p),Ya);
        end
        pi_o = max(0,min(cap(:),pi_o));

        try
            [b,s2] = drift_and_sigma2(Ya,pi_o,p);
            if isscalar(b),  b  = repmat(b,size(Ya));  end
            if isscalar(s2), s2 = repmat(s2,size(Ya)); end
            if numel(b) ~= numel(Ya), error('not vectorized'); end
        catch
            b = zeros(numel(Ya),1);
            s2 = zeros(numel(Ya),1);
            for q = 1:numel(Ya)
                [b(q),s2(q)] = drift_and_sigma2(Ya(q),pi_o(q),p);
            end
        end
        b = b(:);
        s2 = max(0,s2(:));

        Yn = min(Ya+b*dt+sqrt(s2).*sdt.*randn(size(Ya)),ystar);

        % Brownian-bridge correction for crossing the intervention boundary.
        U = rand(numel(Ya),1);
        cr = (Yn <= p.underline_y);
        cand = ~cr & (Ya>p.underline_y) & (Yn>p.underline_y) & (s2>1e-14);
        if any(cand)
            jj = find(cand);
            ex = -2*(Ya(jj)-p.underline_y).*(Yn(jj)-p.underline_y) ./ ...
                (s2(jj)*dt);
            cr(jj) = U(jj) < exp(min(0,ex));
        end

        if any(cr)
            jc = ia(cr);
            nint(jc) = nint(jc)+1;

            if r.active_issue
                iss(jc) = iss(jc)+exp(-rho_L*t)*r.xi_post;
                Yn(cr) = r.y_post;
            else
                alive(jc) = false;
                Yn(cr) = p.underline_y;
            end
        end

        Y(ia) = Yn;
    end

    mc = struct();
    mc.E_N = mean(nint);
    mc.E_issue = mean(iss);
end

function [y_post, xi_post, active] = boundary_post_target(y,v,p)

    y_lo = p.underline_y;
    om = 1-p.kappa_p;

    obj = v - (y-y_lo)/om ...
        - (p.kappa*p.kappa_p/om)*(y_lo-1);
    obj(y<=y_lo) = -Inf;

    [Hval,j] = max(obj);
    active = Hval > (y_lo-1)+1e-8;

    if active
        y_post = y(j);
        xi_post = max( ...
            (y_post-(1-p.kappa)*y_lo-p.kappa)/om, 0);
    else
        y_post = NaN;
        xi_post = NaN;
    end
end
