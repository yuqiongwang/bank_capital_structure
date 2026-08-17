function pim = compute_pi_bar(y, params)
    if y <= 1
        pim = 0;
        return;
    end
    t1 = (1/params.a1) * (1 - 1/y);   % solvency cap
    t2 = (1/params.a3) * (1 - params.a2/y);   % liquidity (LCR) cap
    % constraint regimes
    % both corresponds to the original case
    % none means no regulation at all
    if isfield(params, 'regime') && ~isempty(params.regime)
        regime = params.regime;
    else
        regime = 'both';
    end
    switch regime
        case 'none'        
            if ~isempty(params.pi_cap)
                pim = params.pi_cap;
            else
                pim = 10;   %put a large cap for discretization
            end
        case 'solvency'    
            pim = t1;
        case 'lcr'        
            pim = t2;
        otherwise         
            pim = min(t1, t2);
    end
    pim = max(0, pim);
    if ~isempty(params.pi_cap)
        pim = min(pim, params.pi_cap);
    end
end