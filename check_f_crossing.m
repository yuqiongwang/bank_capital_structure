function [has_cross, fmin, fmax] = check_f_crossing(y, v, pi_max, params)
    rho_L = params.rho - params.mu_L;
    mu_star = params.r + pi_max .* (params.mu - params.r);
    f = rho_L.* v - params.gamma + (params.mu_L - mu_star).* y;
    fmin = min(f);
    fmax = max(f);
    has_cross = (fmin <= 0) && (fmax >= 0);
end

