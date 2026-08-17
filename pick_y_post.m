function [y_post, method] = pick_y_post(y, v, v_prime, y_star, kappa_p)
    % v'(y_post) = 1/(1-kappa_p)
    % which is in the interior of the continuation region
    target = 1 / (1 - kappa_p);
    I = find(y <= y_star);
    if numel(I) < 5
        [y_post, method] = pick_y_post_from_H_discrete(y, v, kappa_p);
        return;
    end
    yI = y(I);
    g  = v_prime(I) - target;
    i0 = min(5, numel(yI)-1);
    y_post = NaN;
    method = 'slope';
    for j = i0+1:numel(yI)
        if g(j-1) > 0 && g(j) <= 0
            denom = g(j-1) - g(j);
            if abs(denom) > 1e-14
                alpha = g(j-1) / denom;
                y_post = yI(j-1) + alpha * (yI(j) - yI(j-1));
            else
                y_post = 0.5 * (yI(j-1) + yI(j));
            end
            break;
        end
    end
    if isnan(y_post)
        [y_post, method] = pick_y_post_from_H_discrete(y, v, kappa_p);
    end
end
