function P_surv = estimate_survival_prob(sol, params, y0, T_horizon, n_paths, dt)
    % first hit of the lower barrier.
    y       = sol.y;
    pi_star = sol.pi_star;
    y_star  = sol.y_star;
    y_lo    = params.underline_y;

    n_steps = ceil(T_horizon / dt);
    sqrt_dt = sqrt(dt);

    survived = 0;
    for p = 1:n_paths
        Y = y0;
        hit = false;
        for step = 1:n_steps
            pi_opt = interp1(y, pi_star, Y, 'linear', 'extrap');
            pi_opt = max(0, min(compute_pi_bar(Y, params), pi_opt));
            [b, sigma2] = drift_and_sigma2_physical(Y, pi_opt, params);
            sigmaY = sqrt(max(0, sigma2));
            Y = Y + b * dt + sigmaY * sqrt_dt * randn();
            if Y > y_star          
                Y = y_star;
            end
            if Y <= y_lo          
                hit = true;
                break;
            end
        end
        if ~hit
            survived = survived + 1;
        end
    end
    P_surv = survived / n_paths;
end
