function sol = solve_bank_qvi(params, settings, mode)
    N = settings.N;
    y_lo = params.underline_y;% the intervention threshold           
    y = linspace(y_lo, settings.y_max, N)';
    dy = y(2) - y(1);
    rhoL = params.rho - params.mu_L;
    inv_dy  = 1/dy;
    inv_dy2 = 1/(dy^2);
    pi_max = zeros(N, 1);
    for i = 1:N
        pi_max(i) = compute_pi_bar(y(i), params);
    end
    %initialize
    v = max(y - 1, 0);
    if strcmp(mode, 'HJB')
        v(1) = y_lo - 1;% where we liquidate, if we don't choose
        % to recapitalize
    end
    pi_star = zeros(N, 1);
    pi_fracs = linspace(0, 1, settings.pi_grid_size);
    M = length(pi_fracs);
    for iter = 1:settings.max_iter
        v_old = v;
        v(N) = 2*v(N-1) - v(N-2);  
        if strcmp(mode, 'HJB')
            v(1) = y_lo - 1;
        end

        %find best \pi
        for i = 2:N-1
            yi = y(i);
            pim = pi_max(i);

            best_val = -1e100;
            best_pi = 0;

            for k = 1:M
                pi = pim * pi_fracs(k);
                [b, s2] = drift_and_sigma2(yi, pi, params);
                A = 0.5 * s2 * inv_dy2;

                bp = max(b, 0);
                bm = max(-b, 0);
                A_im1 = A + bm * inv_dy;
                A_ip1 = A + bp * inv_dy;
                A_i   = -(A_im1 + A_ip1) - rhoL;

                Lv = A_im1 * v(i-1) + A_i * v(i) + A_ip1 * v(i+1);

                if Lv > best_val
                    best_val = Lv;
                    best_pi = pi;
                end
            end

            pi_star(i) = best_pi;
        end
        pi_star(1) = pi_star(2);
        pi_star(N) = pi_star(N-1);

        %solve for v
        n_int = N - 2;
        a_sub = zeros(n_int-1, 1);
        b_diag = zeros(n_int, 1);
        c_sup = zeros(n_int-1, 1);
        rhs = zeros(n_int, 1);

        if strcmp(mode, 'HJB')
            v0 = y_lo - 1;
        else
            v0 = v(1);
        end
        vN = v(N);

        for idx = 1:n_int
            i = idx + 1;
            yi = y(i);
            pi = pi_star(i);

            [b, s2] = drift_and_sigma2(yi, pi, params);
            A = 0.5 * s2 * inv_dy2;

            bp = max(b, 0);
            bm = max(-b, 0);
            A_im1 = A + bm * inv_dy;
            A_ip1 = A + bp * inv_dy;
            A_i   = -(A_im1 + A_ip1) - rhoL;

            b_diag(idx) = 1 - settings.dt * A_i;
            rhs(idx) = v(i);

            if idx > 1
                a_sub(idx-1) = -settings.dt * A_im1;
            else
                rhs(idx) = rhs(idx) - (-settings.dt * A_im1) * v0;
            end

            if idx < n_int
                c_sup(idx) = -settings.dt * A_ip1;
            else
                rhs(idx) = rhs(idx) - (-settings.dt * A_ip1) * vN;
            end
        end

        v_int = thomas_solve(a_sub, b_diag, c_sup, rhs);

        v_new = v;
        v_new(2:N-1) = v_int;
        v_new(N) = v_new(N-1) + dy;

       
        v_new = max(v_new, y - 1);

        if strcmp(mode, 'QVI')
            Hv = compute_H_operator(y, v_new, params);
            v_new = max(v_new, Hv);
            v_new(1) = max(y_lo - 1, Hv(1));
            % here can choose to recapitalize if it is better
        else
            v_new(1) = y_lo - 1;
        end

        %enforce v' >= 1
        v_new = enforce_slope_floor(v_new, dy);
        if strcmp(mode, 'HJB')
            v_new(1) = y_lo - 1;
        end

        v = v_new;

        err = max(abs(v - v_old));
        if err < settings.tol
            break;
        end
    end

    % near y=1
    v_prime  = zeros(N, 1);
    v_second = zeros(N, 1);

    if N >= 4
        v_prime(1) = (-3*v(1) + 4*v(2) - v(3)) / (2*dy);
        v_prime(N) = ( 3*v(N) - 4*v(N-1) + v(N-2)) / (2*dy);

        %centered fd
        v_prime(2:N-1) = (v(3:N) - v(1:N-2)) / (2*dy);

        v_second(1) = (2*v(1) - 5*v(2) + 4*v(3) - v(4)) / (dy^2);
        v_second(N) = (2*v(N) - 5*v(N-1) + 4*v(N-2) - v(N-3)) / (dy^2);
        v_second(2:N-1) = (v(3:N) - 2*v(2:N-1) + v(1:N-2)) / (dy^2);
    else
        v_prime(2:N-1) = (v(3:N) - v(1:N-2)) / (2*dy);
        v_prime(1) = (v(2) - v(1)) / dy;
        v_prime(N) = (v(N) - v(N-1)) / dy;
        v_second(2:N-1) = (v(3:N) - 2*v(2:N-1) + v(1:N-2)) / (dy^2);
        v_second(1) = v_second(2);
        v_second(N) = v_second(N-1);
    end

    %compute y^*
    if strcmp(mode, 'QVI')
        y_star = compute_y_star_analytical(y, v, pi_max, params, pi_star);
    else
        y_star = compute_y_star_from_slope(y, v);
    end

    %compute Hv
    if strcmp(mode, 'QVI')
        Hv = compute_H_operator(y, v, params);
    else
        Hv = -inf(N, 1);
    end

    sol.y = y;
    sol.v = v;
    sol.v_prime = v_prime;
    sol.v_second = v_second;
    sol.pi_star = pi_star;
    sol.pi_max = pi_max;
    sol.y_star = y_star;
    sol.Hv = Hv;
    sol.mode = mode;
    sol.iterations = iter;
end
