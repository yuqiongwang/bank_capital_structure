function x = thomas_solve(a, b, c, d)
%thomas solver
    n = length(b);
    cp = zeros(n-1, 1);
    dp = zeros(n, 1);
    cp(1) = c(1) / b(1);
    dp(1) = d(1) / b(1);
    for i = 2:n-1
        denom = b(i) - a(i-1) * cp(i-1);
        cp(i) = c(i) / denom;
        dp(i) = (d(i) - a(i-1) * dp(i-1)) / denom;
    end
    dp(n) = (d(n) - a(n-1) * dp(n-1)) / (b(n) - a(n-1) * cp(n-1));
    x = zeros(n, 1);
    x(n) = dp(n);
    for i = n-1:-1:1
        x(i) = dp(i) - cp(i) * x(i+1);
    end
end
