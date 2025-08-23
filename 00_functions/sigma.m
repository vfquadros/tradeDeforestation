function Sigma = sigma(params,S)
        Sigma = sqrt(max(0,params.Sbar^2 - params.elasticity.*(params.Sbar - S).^2));
        Sigma = max(0,sqrt(max(0,params.Sbar^2 - params.elasticity.*(params.Sbar - S).^2)));

end 