function check_params_valid(params)
    r = params.r; mu = params.mu; mu_L = params.mu_L;
    rho = params.rho; gamma = params.gamma;
    a1 = params.a1; a2 = params.a2; a3 = params.a3;
    mu_r_plus = max(0, mu - r);
    a_bar = max(a1, a3);
    rho_L = rho - mu_L;
    threshold1 = max(mu_L, r + mu_r_plus / a_bar);
    if rho > threshold1
        fprintf('condition 1: ok\n');
    else
        fprintf('condition 1: not ok\n');
    end
    if abs(a3 - a1) < 1e-10, y_hat = 1;
    else, y_hat = (a3 - a1*a2) / (a3 - a1); end
    term_A = r + mu_r_plus/a3 - rho;
    A = term_A * y_hat - a2 * mu_r_plus / a3;
    term_B = r + mu_r_plus/a1 - rho;
    B = max(0, term_B)*y_hat - max(0, -term_B) - mu_r_plus/a1;
    threshold2 = max(A, B) + gamma;
    if -rho_L < threshold2
        fprintf('condition 2: ok\n');
    else
        fprintf('Condition 2: not ok\n');
    end
end
