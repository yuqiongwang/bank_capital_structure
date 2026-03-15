function metrics = compute_metrics(sol_qvi, sol_hjb)
    y = sol_qvi.y;
    [~, idx_12] = min(abs(y - 1.2));
    metrics.y_star_qvi = sol_qvi.y_star;
    metrics.y_star_hjb = sol_hjb.y_star;
    %
    metrics.v_at_1_qvi  = sol_qvi.v(1);
    metrics.v_at_1_hjb  = sol_hjb.v(1);
    %
    metrics.v_at_12_qvi = sol_qvi.v(idx_12);
    metrics.v_at_12_hjb = sol_hjb.v(idx_12);
    %
    metrics.Delta_v_at_1  = metrics.v_at_1_qvi  - metrics.v_at_1_hjb;
    metrics.Delta_v_at_12 = metrics.v_at_12_qvi - metrics.v_at_12_hjb;
    metrics.Delta_v_max   = max(sol_qvi.v - sol_hjb.v);
end