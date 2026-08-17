function [ok, info] = check_params_valid(params)
% check the well-posedness of the parameters and return 
% in each case if it is ok
% this incluse both the wellposeness, exclude the immediate liquidation 
% case, and ask r>r_L
% purely for numerical experiments purposes, the latter two are not
% strictly necessary
    r = params.r; mu = params.mu; mu_L = params.mu_L;
    rho = params.rho; gamma = params.gamma;
    a1 = params.a1; a2 = params.a2; a3 = params.a3;

    mu_r_plus = max(0, mu - r);
    a_bar = max(a1, a3);
    rho_L = rho - mu_L;

    if isfield(params, 'pi_cap') && ~isempty(params.pi_cap)
        pi_sup = min(1/a_bar, params.pi_cap);
    else
        pi_sup = 1/a_bar;
    end
    sup_drift = r + mu_r_plus * pi_sup;

    thr1    = max(mu_L, sup_drift);
    c1      = rho > thr1;
    margin1 = rho - thr1;

    if abs(a3 - a1) < 1e-10, y_hat = 1;
    else, y_hat = (a3 - a1*a2) / (a3 - a1); end
    term_A = r + mu_r_plus/a3 - rho;
    A = term_A * y_hat - a2 * mu_r_plus / a3;
    term_B = r + mu_r_plus/a1 - rho;
    B = max(0, term_B)*y_hat - max(0, -term_B) - mu_r_plus/a1;
    thr2    = max(A, B) + gamma;
    c2      = (-rho_L) < thr2;
    margin2 = thr2 - (-rho_L);

    % r>r_L is not strictly necessary for the wellposeness, 
    % but we put here as a simulation check
    r_L  = mu_L - gamma;
    c_rL = r > r_L;

    ok = c1 && c2 && c_rL;

    info = struct('c1',c1,'c2',c2,'c_rL',c_rL, ...
                  'margin1',margin1,'margin2',margin2, ...
                  'margin',min([margin1, margin2]), ...
                  'thr1',thr1,'thr2',thr2,'pi_sup',pi_sup, ...
                  'r_L',r_L,'a_bar',a_bar);

    if nargout == 0
        fprintf('C1 (wellposeness)   : %s (margin=%+.4f)\n', tf2str(c1), margin1);
        fprintf('C2 (no liquidation) : %s (margin=%+.4f)\n', tf2str(c2), margin2);
        fprintf('r > r_L             : %s (r=%.4f, r_L=%.4f)\n', tf2str(c_rL), r, r_L);
    end
end

function s = tf2str(tf)
    if tf, s = 'ok'; else, s = 'FAIL'; end
end
