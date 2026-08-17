%% Dense-trigger grid for coordinated and stackelberg 
% we compare in this case what is the difference in \underline y^B
% when the regulator chooses the threshold, versus the bank chooses the
% threshold
clear; clc; close all;
%% 
use_parallel    = true;
chunk_size      = 40;    
force_restart   = false;

%% baseline parameters
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
base_params.kappa_p = 0.02;
base_params.regime  = 'both';
base_params.a1 = 0.045; base_params.a2 = 0.05; base_params.a3 = 0.30;
base_params.underline_y = 1.02;

%% scenarios 
scenario_name   = ["tight_cap"; "loose_cap"];
scenario_pi_cap = [1.0;          10];

%% mesh and domain sizing
h_master     = 0.005;
y_max_safety = 1.15;  
y_max_pad    = 0.50;  
y_max_hard   = 8.00;

%% grids
a1_grid = [0.045, 0.060, 0.080, 0.100, 0.120, 0.150];
a2_grid = [0.050, 0.100, 0.150];
a3_grid = [0.150, 0.300, 0.500];

uly_step_by_scenario = [0.005; 0.010];   % tight; loose
uly_grids = cell(numel(scenario_name),1);
for is = 1:numel(scenario_name)
    st = uly_step_by_scenario(is);
    uly_grids{is} = (round(1.010/st) : round(1.180/st)) * st;
    check_on_lattice(uly_grids{is}, h_master);
end

T_horizon = 5.0;
y0        = 1.20;
y_ref     = max(cellfun(@max, uly_grids));

check_on_lattice([y0, y_ref, y_max_hard], h_master);

%% solver settings
settings = struct();
settings.dt           = 0.05;
settings.max_iter     = 3000;
settings.tol          = 1e-9;
settings.pi_grid_size = 151;
settings.verbose      = false;
settings.compute_analytical_ystar = false;

%% Monte Carlo
n_paths_screen        = 5000;
dt_mc                 = 0.05;
rng_seed              = 12345;
alpha_ci              = 0.05;
use_bridge_correction = true;
track_interventions   = true;

%%
n_scen  = numel(scenario_name);
n_total = 0;
for is = 1:n_scen
    n_total = n_total + numel(a1_grid)*numel(a2_grid)*numel(a3_grid)*numel(uly_grids{is});
end

R      = init_storage(n_total);
rowkey = strings(n_total, 1);
kmap   = containers.Map('KeyType','char','ValueType','double');
Wkey = {}; Wp = {}; Ws = {};
idx = 0;

for is = 1:n_scen
  uly_grid = uly_grids{is};
  for i1 = 1:numel(a1_grid)
    for i2 = 1:numel(a2_grid)
      for i3 = 1:numel(a3_grid)
        for iu = 1:numel(uly_grid)

          idx = idx + 1;

          p = base_params;
          p.pi_cap      = scenario_pi_cap(is);
          p.a1          = a1_grid(i1);
          p.a2          = a2_grid(i2);
          p.a3          = a3_grid(i3);
          p.underline_y = uly_grid(iu);

          R.row_id(idx)       = idx;
          R.scenario(idx)     = scenario_name(is);
          R.pi_cap(idx)       = p.pi_cap;
          R.kappa(idx)        = p.kappa;
          R.kappa_p(idx)      = p.kappa_p;
          R.a1(idx)           = p.a1;
          R.a2(idx)           = p.a2;
          R.a3(idx)           = p.a3;
          R.underline_y(idx)  = p.underline_y;
          R.stress_level(idx) = p.underline_y;
          R.y_ref(idx)        = y_ref;

          [ok, info] = assumption31_holds(p);
          R.assumption_margin(idx) = info.margin;
          R.y_star_bound(idx)      = info.y_star_bound;

          if ~ok
              R.message(idx) = "Assumption 3.1 fails";
              R.done(idx) = true;
              continue;
          end

          y_max_target = max(y_max_safety*info.y_star_bound, ...
                             p.underline_y + y_max_pad);

          if y_max_target > y_max_hard
              R.message(idx) = "y* bound exceeds y_max_hard";
              R.done(idx) = true;
              R.domain_warning(idx) = true;
              continue;
          end

          y_max = p.underline_y + ...
              ceil((y_max_target - p.underline_y)/h_master) * h_master;

          s = settings;
          s.y_max = y_max;
          s.N     = round((y_max - p.underline_y)/h_master) + 1;

          R.y_max(idx) = y_max;
          R.N(idx)     = s.N;
          R.h(idx)     = h_master;

          inert = lcr_is_dominated(p, y_max);
          R.lcr_inert(idx) = inert;

          if inert
              key = sprintf('%s|%.6f|%.6f|%.6f', ...
                  scenario_name(is), p.pi_cap, p.a1, p.underline_y);
          else
              key = sprintf('%s|%.6f|%.6f|%.6f|%.6f|%.6f', ...
                  scenario_name(is), p.pi_cap, p.a1, p.a2, p.a3, p.underline_y);
          end

          rowkey(idx) = string(key);

          if ~isKey(kmap, key)
              kmap(key) = numel(Wkey) + 1;
              Wkey{end+1} = key;   %#ok<SAGROW>
              Wp{end+1}   = p;     %#ok<SAGROW>
              Ws{end+1}   = s;     %#ok<SAGROW>
          end
        end
      end
    end
  end
end

nW = numel(Wkey);
Nw = cellfun(@(x) x.N, Ws);

for is = 1:n_scen
    fprintf('%-10s triggers : %d (step %.4f)\n', ...
        scenario_name(is), numel(uly_grids{is}), uly_step_by_scenario(is));
end
fprintf('Rows enumerated : %d\n', n_total);
fprintf('Rows skipped    : %d\n', nnz(R.done));
fprintf('Unique solves   : %d\n', nW);
fprintf('N per solve     : min %d, median %d, max %d\n', ...
    min(Nw), round(median(Nw)), max(Nw));
for is = 1:n_scen
    c = sum(startsWith(string(Wkey), scenario_name(is)));
    fprintf('  %-10s : %d solves\n', scenario_name(is), c);
end
fprintf('\n');

cfg_hash = sum(double([Wkey{:}]));   
res    = cell(nW,1);
okv    = false(nW,1);
msg    = strings(nW,1);
done_w = false(nW,1);

has_pct = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));

if use_parallel && has_pct
    if isempty(gcp('nocreate'))
        try
            parpool;
        catch
            warning('Could not start a pool; running serially.');
            use_parallel = false;
        end
    end
elseif use_parallel
    fprintf('Parallel Computing Toolbox not available; running serially.\n\n');
    use_parallel = false;
end

todo     = find(~done_w);
n_chunks = max(ceil(numel(todo)/chunk_size), 0);
t_start  = tic;
n_start  = nnz(done_w);

for ic = 1:n_chunks
    lo = (ic-1)*chunk_size + 1;
    hi = min(ic*chunk_size, numel(todo));
    w_list = todo(lo:hi);

    cWp = Wp(w_list);
    cWs = Ws(w_list);
    m   = numel(w_list);
    cres = cell(m,1); cok = false(m,1); cmsg = strings(m,1);

    if use_parallel
        parfor q = 1:m
            try
                cres{q} = solve_and_simulate(cWp{q}, cWs{q}, y0, y_ref, ...
                    T_horizon, n_paths_screen, dt_mc, rng_seed, ...
                    use_bridge_correction, track_interventions);
                cok(q) = true; cmsg(q) = "ok";
            catch ME
                cmsg(q) = string(ME.message);
            end
        end
    else
        for q = 1:m
            try
                cres{q} = solve_and_simulate(cWp{q}, cWs{q}, y0, y_ref, ...
                    T_horizon, n_paths_screen, dt_mc, rng_seed, ...
                    use_bridge_correction, track_interventions);
                cok(q) = true; cmsg(q) = "ok";
            catch ME
                cmsg(q) = string(ME.message);
            end
        end
    end

    res(w_list)    = cres;
    okv(w_list)    = cok;
    msg(w_list)    = cmsg;
    done_w(w_list) = true;

    el    = toc(t_start);
    ndone = nnz(done_w) - n_start;
    rate  = el / max(ndone, 1);
    left  = (nW - nnz(done_w)) * rate;
    fprintf(['chunk %3d/%3d | %4d/%4d solves | %6.1f min elapsed | ', ...
             '%.2f s/solve | ~%.0f min left\n'], ...
        ic, n_chunks, nnz(done_w), nW, el/60, rate, left/60);
end

nfail = nnz(done_w & ~okv);
if nfail > 0
    warning('%d solves failed; see the message column.', nfail);
end

for i = 1:n_total
    if R.done(i) || rowkey(i) == "", continue; end
    w = kmap(char(rowkey(i)));
    if okv(w)
        R = write_result(R, i, res{w}, alpha_ci, n_paths_screen);
        R.reused(i) = true;
    else
        R.message(i) = msg(w);
        R.done(i)    = true;
    end
end

T = storage_to_table(R);
T = T(T.done | T.valid, :);


%% local functions
function check_on_lattice(vals, h)
    bad = vals(abs(round(vals/h) - vals/h) > 1e-9);
    if ~isempty(bad)
        error('Not multiples of h=%.6f: %s', h, mat2str(bad));
    end
end

function tf = lcr_is_dominated(p, y_max)
    yy  = linspace(p.underline_y, y_max, 400);
    lcr = (1./p.a3) .* (1 - p.a2 ./ yy);
    sol = (1./p.a1) .* (1 - 1 ./ yy);
    tf  = all(lcr >= min(sol, p.pi_cap) - 1e-12);
end

function res = solve_and_simulate(p, s, y0, y_ref, T_horizon, ...
        n_paths, dt, seed, use_bridge, track)

    evalc('sol = solve_bank_qvi(p, s, ''QVI'');');

    [d0, j0] = min(abs(sol.y - y0));
    if d0 > 1e-9
        error('y0=%.4f is not a grid node (nearest %.6f).', y0, sol.y(j0));
    end

    res = struct();
    res.v_y0   = sol.v(j0);
    res.y_star = sol.y_star;
    res.y_max  = s.y_max;
    res.h      = sol.y(2) - sol.y(1);

    if isfield(sol, 'n_iter'), res.n_iter = sol.n_iter; else, res.n_iter = NaN; end

    [y_post, xi_post, active_issue] = get_boundary_post_target(sol.y, sol.v, p);
    res.y_post       = y_post;
    res.xi_post      = xi_post;
    res.active_issue = active_issue;

    res.domain_warning = ~isfinite(sol.y_star) || (sol.y_star > s.y_max - 20*res.h);

    jump = struct('active', active_issue, 'y_post', y_post, 'xi_post', xi_post);

    M = simulate_metrics(sol, p, jump, y0, y_ref, T_horizon, n_paths, dt, ...
        seed, use_bridge, track);

    res.P_own = M.P_own;  res.k_own   = M.k_own;
    res.P_ref = M.P_ref;  res.k_ref   = M.k_ref;
    res.E_N   = M.E_N;    res.E_issue = M.E_issue;
end

function M = simulate_metrics(sol, p, jump, y0, y_ref, T_horizon, n_paths, dt, ...
        seed, use_bridge, track)

    rng(seed, 'twister');

    yg = sol.y(:); pig = sol.pi_star(:); ystar = sol.y_star;
    rho_L = p.rho - p.mu_L;

    n_steps = ceil(T_horizon/dt);
    sdt = sqrt(dt);

    Y        = y0*ones(n_paths,1);
    hit_own  = false(n_paths,1);
    hit_ref  = false(n_paths,1);
    nint     = zeros(n_paths,1);
    disc_iss = zeros(n_paths,1);
    alive    = true(n_paths,1);

    ref_meaningful = (p.underline_y <= y_ref + 1e-12);

    for step = 1:n_steps
        ia = find(alive);
        if isempty(ia), break; end
        t = step*dt;

        Ya = Y(ia);

        pi_o = interp1(yg, pig, Ya, 'linear', 'extrap');
        try
            cap = compute_pi_bar(Ya, p);
            if isscalar(cap), cap = repmat(cap, size(Ya)); end
            if numel(cap) ~= numel(Ya), error('not vectorized'); end
        catch
            cap = arrayfun(@(yy) compute_pi_bar(yy, p), Ya);
        end
        pi_o = max(0, min(cap(:), pi_o));

        try
            [b, s2] = drift_and_sigma2(Ya, pi_o, p);
            if isscalar(b),  b  = repmat(b,  size(Ya)); end
            if isscalar(s2), s2 = repmat(s2, size(Ya)); end
            if numel(b) ~= numel(Ya), error('not vectorized'); end
        catch
            b = zeros(numel(Ya),1); s2 = zeros(numel(Ya),1);
            for q = 1:numel(Ya)
                [b(q), s2(q)] = drift_and_sigma2(Ya(q), pi_o(q), p);
            end
        end
        b = b(:); s2 = max(0, s2(:));

        Yn = Ya + b*dt + sqrt(s2).*sdt.*randn(size(Ya));
        Yn = min(Yn, ystar);
        U = rand(numel(Ya), 1);

        if ref_meaningful
            cr = cross_indicator(Ya, Yn, y_ref, s2, dt, use_bridge, U);
            hit_ref(ia(cr)) = true;
        end

        co = cross_indicator(Ya, Yn, p.underline_y, s2, dt, use_bridge, U);
        hit_own(ia(co)) = true;

        if any(co)
            jc = ia(co);
            nint(jc) = nint(jc) + 1;

            if track && jump.active
                disc_iss(jc) = disc_iss(jc) + exp(-rho_L*t)*jump.xi_post;
                Yn(co) = jump.y_post;
            else
                alive(jc) = false;
                Yn(co) = p.underline_y;
            end
        end

        Y(ia) = Yn;
    end

    M = struct();
    M.k_own = nnz(~hit_own);
    M.P_own = M.k_own / n_paths;

    if ref_meaningful
        M.k_ref = nnz(~hit_ref);
        M.P_ref = M.k_ref / n_paths;
    else
        M.k_ref = NaN; M.P_ref = NaN;
    end

    M.E_N     = mean(nint);
    M.E_issue = mean(disc_iss);
end

function cr = cross_indicator(Ya, Yn, level, s2, dt, use_bridge, U)
    cr = (Yn <= level);
    if use_bridge
        cand = ~cr & (Ya > level) & (Yn > level) & (s2 > 1e-14);
        if any(cand)
            j = find(cand);
            ex = -2*(Ya(j)-level).*(Yn(j)-level) ./ (s2(j)*dt);
            cr(j) = U(j) < exp(min(0, ex));
        end
    end
end

function [y_post, xi_post, active_issue] = get_boundary_post_target(y, v, p)
    y_lo = p.underline_y;
    om = 1 - p.kappa_p;
    obj = v - (y - y_lo)/om - (p.kappa*p.kappa_p/om)*(y_lo - 1);
    obj(y <= y_lo) = -Inf;
    [Hval, j] = max(obj);
    liq = y_lo - 1;
    active_issue = Hval > liq + 1e-8;
    if active_issue
        y_post  = y(j);
        xi_post = max((y_post - (1-p.kappa)*y_lo - p.kappa)/om, 0);
    else
        y_post = NaN; xi_post = NaN;
    end
end

function [lo, hi] = wilson_ci(k, n, alpha)
    if n <= 0 || ~isfinite(k), lo = NaN; hi = NaN; return; end
    z = -sqrt(2)*erfcinv(2*(1-alpha/2));
    ph = k/n; den = 1 + z^2/n;
    c  = (ph + z^2/(2*n))/den;
    hw = z*sqrt((ph*(1-ph) + z^2/(4*n))/n)/den;
    lo = max(0, c-hw); hi = min(1, c+hw);
end

function R = init_storage(n)
    z = @() NaN(n,1);
    R = struct( ...
        'row_id',z(),'scenario',strings(n,1),'pi_cap',z(), ...
        'kappa',z(),'kappa_p',z(), ...
        'a1',z(),'a2',z(),'a3',z(),'underline_y',z(),'stress_level',z(),'y_ref',z(), ...
        'y_star',z(),'y_post',z(),'xi_post',z(),'active_issue',false(n,1), ...
        'v_y0',z(),'y_max',z(),'N',z(),'h',z(),'n_iter',z(), ...
        'P_no_intervention',z(),'P_no_intervention_lo',z(),'P_no_intervention_hi',z(), ...
        'n_success',z(), ...
        'P_ref',z(),'P_ref_lo',z(),'P_ref_hi',z(), ...
        'E_N',z(),'E_disc_issuance',z(), ...
        'assumption_margin',z(),'y_star_bound',z(), ...
        'lcr_inert',false(n,1),'reused',false(n,1),'domain_warning',false(n,1), ...
        'valid',false(n,1),'done',false(n,1),'message',strings(n,1));
end

function R = write_result(R, idx, res, alpha, npaths)
    R.v_y0(idx)    = res.v_y0;
    R.y_star(idx)  = res.y_star;
    R.y_post(idx)  = res.y_post;
    R.xi_post(idx) = res.xi_post;
    R.n_iter(idx)  = res.n_iter;
    R.active_issue(idx)   = res.active_issue;
    R.domain_warning(idx) = res.domain_warning;

    R.P_no_intervention(idx) = res.P_own;
    R.n_success(idx) = res.k_own;
    [lo, hi] = wilson_ci(res.k_own, npaths, alpha);
    R.P_no_intervention_lo(idx) = lo;
    R.P_no_intervention_hi(idx) = hi;

    R.P_ref(idx) = res.P_ref;
    [lo, hi] = wilson_ci(res.k_ref, npaths, alpha);
    R.P_ref_lo(idx) = lo;
    R.P_ref_hi(idx) = hi;

    R.E_N(idx) = res.E_N;
    R.E_disc_issuance(idx) = res.E_issue;

    R.valid(idx) = true;
    R.done(idx)  = true;
    R.message(idx) = "ok";
end

function T = storage_to_table(R)
    f = fieldnames(R);
    args = cell(1, numel(f));
    for k = 1:numel(f), args{k} = R.(f{k}); end
    T = table(args{:}, 'VariableNames', f');
end
