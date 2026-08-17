%% robustness of the frontier across horizon T and initial health y0
% this is not disclosed numerically in the paper, we just check if 
% the result is robust
clear; clc; close all;
%% parameters (baseline, r > r_L)
base_params.r = 0.02; base_params.mu = 0.04; base_params.mu_L = 0.03;
base_params.sigma = 0.08; base_params.sigma_L = 0.03; base_params.c = 0.20;
base_params.gamma = 0.02; base_params.rho = 0.12;
base_params.kappa = 0.01; base_params.kappa_p = 0.02;
base_params.pi_cap = 1;
base_params.a1 = 0.045; base_params.a2 = 0.05; base_params.a3 = 0.30;
base_params.underline_y = 1.02;

a1_grid = [0.045, 0.06, 0.08, 0.10, 0.12];
a2_grid = [0.05, 0.10, 0.18];
a3_grid = [0.15, 0.30, 0.50];

T_list  = [1.0, 5.0, 10.0];       
y0_list = [1.10, 1.20, 1.30];% y=1.2 is tbe baseline
settings = struct('N',101,'y_max',2.5,'dt',0.05,'max_iter',1500,'tol',1e-6,'pi_grid_size',151);
n_paths = 5000; dt = 0.05;
mc_seed = 20240517; 

n1=numel(a1_grid); n2=numel(a2_grid); n3=numel(a3_grid); n_triples=n1*n2*n3;
A1=zeros(n_triples,1); A2=A1; A3=A1; YS=A1;
sols = cell(n_triples,1);
t=0;
for i1=1:n1, for i2=1:n2, for i3=1:n3
    t=t+1; p=base_params; p.a1=a1_grid(i1); p.a2=a2_grid(i2); p.a3=a3_grid(i3);
    A1(t)=p.a1; A2(t)=p.a2; A3(t)=p.a3;
    [ok_feas, ~] = check_params_valid(p);  
    if ~ok_feas, sols{t}=[]; YS(t)=NaN; continue; end
    try
        s = solve_bank_qvi(p, settings, 'QVI'); sols{t}=s; YS(t)=s.y_star;
    catch
        sols{t}=[]; YS(t)=NaN;
    end
end, end, end

%% loop over (T, y0)
fprintf('\n%-6s %-6s | %-8s | frontier size | share with (a2,a3)=baseline\n', 'T','y0','corr');
summary = {};
for it = 1:numel(T_list)
 for iy = 1:numel(y0_list)
   T = T_list(it); y0 = y0_list(iy);
   v_y0 = nan(n_triples,1); P = nan(n_triples,1);
   for t = 1:n_triples
       s = sols{t};
       if isempty(s) || y0 <= base_params.underline_y, continue; end
       p = base_params; p.a1=A1(t); p.a2=A2(t); p.a3=A3(t);
       v_y0(t) = interp1(s.y, s.v, y0, 'linear', 'extrap');
       rng(mc_seed);   % CRN across triples within each (T,y0) cell
       P(t)    = estimate_survival_prob(s, p, y0, T, n_paths, dt);
   end
   valid = ~isnan(v_y0) & ~isnan(P);
   % Pareto frontier
   is_par=false(n_triples,1); vi=find(valid);
   for ii=1:numel(vi)
       i=vi(ii); dom=false;
       for jj=1:numel(vi)
           j=vi(jj);
           if j~=i && v_y0(j)>=v_y0(i) && P(j)>=P(i) && (v_y0(j)>v_y0(i)||P(j)>P(i))
               dom=true; break;
           end
       end
       is_par(i)=~dom;
   end
   fidx = find(is_par);
   % does the frontier keep (a2,a3) at baseline (i.e. a1 is the active lever)?
   base_a2a3 = (A2(fidx)==base_params.a2) & (A3(fidx)==base_params.a3);
   share = mean(base_a2a3);
   % rank correlation between a1 and P along the frontier (Spearman-like sign)
   if numel(fidx) >= 2
       % Spearman corr; fall back to NaN if Statistics Toolbox absent
       try
           cc = corr(A1(fidx), P(fidx), 'type','Spearman');
       catch
           cc = NaN;
       end
   else
       cc = NaN;
   end
   fprintf('%-6.1f %-6.2f | %-8.3f | %-13d | %.2f\n', T, y0, cc, numel(fidx), share);
   summary(end+1,:) = {T, y0, numel(fidx), share, cc}; %#ok<SAGROW>
 end
end
S = cell2table(summary, 'VariableNames', {'T','y0','frontier_size','a1_only_share','corr_a1_P'});
fprintf('\nRobustness summary (a1_only_share near 1 => a1 still drives the frontier):\n');
disp(S);
