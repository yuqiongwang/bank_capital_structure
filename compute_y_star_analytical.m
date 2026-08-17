function y_star = compute_y_star_analytical(y, v, pi_max, params, pi_star)
    % y* solves  rho_L v(y) = gamma - (mu_L - mu_star(y)) y,
    rho_L = params.rho - params.mu_L;
    N = length(y);

    f = zeros(N, 1);
    for i = 1:N
        if params.mu >= params.r
            mu_opt = params.r + pi_max(i) * (params.mu - params.r);
        else
            mu_opt = params.r;            % pi*_myopic = 0
        end
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
       % fprintf('no root found');
    end
    if y_star <= params.underline_y
        fprintf('y^* <= underline_y');
        y_star = params.underline_y + (y(2) - y(1));
    end
end
