function Vu = value_u(params,A,R,F) 

    term1 = params.alpha .* A .* (1 - params.beta^params.Tlag) ./ (1 - params.beta);
    term2 = params.beta^params.Tlag * (params.pi .* (R+params.wage/(1-params.beta)) + (1 - params.pi) .* ...
        (A + params.beta * params.pi .* (R + params.wage / (1 - params.beta))) ./ (1 - params.beta * (1 - params.pi)));
    Vu = term1 - F + term2;
end
