function y_star = compute_y_star_analytical(y, v, pi_max, params, pi_star)
    % used in VI, pi_max is not used here..
    rho_L = params.rho - params.mu_L;
    N = length(y);

    f = zeros(N, 1);
    for i = 1:N
        mu_opt = params.r + pi_star(i) * (params.mu - params.r);
        f(i) = rho_L * v(i) - params.gamma + (params.mu_L - mu_opt) * y(i);
    end

    y_star = NaN;
    for i = 2:N
        if f(i-1) <= 0 && f(i) >= 0
            if f(i) - f(i-1) > 1e-12
                alpha = -f(i-1) / (f(i) - f(i-1));
                y_star = y(i-1) + alpha * (y(i) - y(i-1));
            else
                y_star = 0.5 * (y(i-1) + y(i));
            end
            break;
        end
    end

    if isnan(y_star)
        [~, idx_min] = min(abs(f));
        y_star = y(idx_min);
        fprintf('no root found');
    end
    if y_star <= 1
        fprintf('y^*<1');
        y_star = 1.01;
    end
end
