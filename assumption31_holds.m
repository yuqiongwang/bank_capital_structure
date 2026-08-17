function [tf, info] = assumption31_holds(p)
% check the well-posedness result Assumption 1
% this has to hold for the analysis of the problem
    a_bar = max(p.a1, p.a3);

    if isfield(p, 'pi_cap') && ~isempty(p.pi_cap)
        pi_inf = min(1 / a_bar, p.pi_cap);
    else
        pi_inf = 1 / a_bar;
    end

    mu_star_sup = p.r + max(p.mu - p.r, 0) * pi_inf;
    threshold   = max(p.mu_L, mu_star_sup);
    margin      = p.rho - threshold;

    tf = margin > 1e-10;

    rho_L = p.rho - p.mu_L;

    if p.rho - mu_star_sup > 1e-12
        y_star_bound = (rho_L + p.gamma) / (p.rho - mu_star_sup);
    else
        y_star_bound = Inf;
    end

    info = struct( ...
        'a_bar',        a_bar, ...
        'pi_inf',       pi_inf, ...
        'mu_star_sup',  mu_star_sup, ...
        'threshold',    threshold, ...
        'margin',       margin, ...
        'rho_L',        rho_L, ...
        'y_star_bound', y_star_bound);
end
