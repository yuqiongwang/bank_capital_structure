function Hv = compute_H_operator(y, v, params)
    kappa = params.kappa;
    kappa_p = params.kappa_p;
    one_minus_kp = 1 - kappa_p;
    N = length(y);
    base = v - y / one_minus_kp;
    suffix = zeros(N, 1);
    suffix(N) = base(N);
    for j = N-1:-1:1
        suffix(j) = max(base(j), suffix(j+1));
    end
    Hv = -1e100 * ones(N, 1);
    for i = 1:N-1
        Hv(i) = suffix(i+1) + ((1-kappa)*y(i) + kappa) / one_minus_kp + kappa * (y(i) - 1);
    end
    Hv(N) = -1e100;
end
