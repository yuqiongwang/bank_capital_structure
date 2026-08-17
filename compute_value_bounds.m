function bounds = compute_value_bounds(params)
    % Upper/lower bounds of v(y) from Prop 3.2
    r     = params.r;
    mu    = params.mu;
    mu_L  = params.mu_L;
    rho   = params.rho;
    gamma = params.gamma;
    rho_L = rho - mu_L;

    z = linspace(params.underline_y, max(params.underline_y + 50, 100), 20001)';
    if mu >= r
        pim = zeros(size(z));
        for i = 1:numel(z)
            pim(i) = compute_pi_bar(z(i), params);
        end
        mu_star = r + pim .* (mu - r);
    else
        mu_star = r * ones(size(z));
    end
    M = max((mu_star - rho) .* z);

    C = max(-rho_L, M + gamma);   

    bounds.rho_L       = rho_L;
    bounds.M           = M;
    bounds.C           = C;
    bounds.upper_const = C / rho_L;   % v(y) <= y + this constant K
end
