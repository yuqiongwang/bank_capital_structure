function y_star = compute_y_star_from_slope(y, v)
    % used in hjb, ask slope hits 1
    dy = y(2)-y(1);
    s  = (v(2:end) - v(1:end-1)) / dy;   % discrete slope
    g  = s - 1;

    y_star = NaN;

    % look for first crossing g>0 to g<=0
    for k = 2:numel(g)
        if g(k-1) > 0 && g(k) <= 0
            denom = g(k-1) - g(k);
            if abs(denom) > 1e-14
                alpha = g(k-1) / denom;              
                y_star = y(k) + alpha*dy;            
            else
                y_star = 0.5*(y(k) + y(k+1));
            end
            break;
        end
    end

    if isnan(y_star)
        fprintf('y*_HJB is probably large');
        % comment: here y_HJB^* large is ok, as we do not use it, and it
        % does not really affect the reported value.
    elseif y_star <= 1
        y_star = 1.01;
    end
end

