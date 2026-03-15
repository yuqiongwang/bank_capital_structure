function [b, s2] = drift_and_sigma2(y, pi, params)
    mu_pi = (1 - pi) * params.r + pi * params.mu;
    b = y * (mu_pi - params.mu_L) + params.gamma;
    s2 = (pi^2) * (params.sigma^2) * (y^2) ...
       + 2* pi * params.c * params.sigma * params.sigma_L * y * (1-y) ...
       + (params.sigma_L^2) * ((1-y)^2);
    s2 = max(0, s2);
end