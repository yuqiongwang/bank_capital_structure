function pim = compute_pi_bar(y, params)
    if y <= 1
        pim = 0;
        return;
    end
    t1 = (1/params.a1) * (1 - 1/y);
    t2 = (1/params.a3) * (1 - params.a2/y);
    pim = max(0, min(t1, t2));
    if ~isempty(params.pi_cap)
        pim = min(pim, params.pi_cap);
    end
end