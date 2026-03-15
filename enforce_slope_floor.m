function v = enforce_slope_floor(v, dy)
    for i = 2:length(v)
        floor_i = v(i-1) + dy;
        if v(i) < floor_i
            v(i) = floor_i;
        end
    end
end