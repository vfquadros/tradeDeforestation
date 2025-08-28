function Profit = profit(params,A,S,k)
        
        gamma = 1 ./ ((params.elasticity-1) .* (1 - params.beta));
        Profit = gamma.*params.delta.*(params.tau.*params.wage./b(params,A,S,k)).^(1-params.elasticity);