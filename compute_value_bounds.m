function bounds = compute_value_bounds(params)
    %upper and lowerbounds of v(y)
    %y - 1 ≤ v(y) ≤ y + (1/ρ_L) * max(-ρ_L, A+γ, B+γ)
    r = params.r;
    mu = params.mu;
    mu_L = params.mu_L;
    rho = params.rho;
    gamma = params.gamma;
    a1 = params.a1;
    a2 = params.a2;
    a3 = params.a3;
    mu_r_plus = max(0, mu - r);
    %
    rho_L = rho - mu_L;
    if abs(a3 - a1) < 1e-10
        y_hat = 1;
    else
        y_hat = (a3 - a1 * a2) / (a3 - a1);
    end
    
    term_A = r + mu_r_plus / a3 - rho;
    A = term_A * y_hat - a2 * mu_r_plus / a3;
    
    term_B = r + mu_r_plus / a1 - rho;
    term_B_plus = max(0, term_B);
    term_B_minus = max(0, -term_B);
    B = term_B_plus * y_hat - term_B_minus - mu_r_plus / a1;
    
    C = max([-rho_L, A + gamma, B + gamma]);
    
    bounds.rho_L = rho_L;
    bounds.y_hat = y_hat;
    bounds.A = A;
    bounds.B = B;
    bounds.C = C;
    bounds.upper_const = C / rho_L;  % v(y) ≤ y + this constant
end