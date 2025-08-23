function Profit = profit(params,A,S,k)
        % k == 1 refers to the profit of plot n conditional to having S
        % other contiguous plots.
        gamma = 1 ./ ((params.elasticity-1) .* (1 - params.beta));
        Profit = gamma.*params.delta.*(params.tau.*params.wage./b(params,A,S,k)).^(1-params.elasticity);