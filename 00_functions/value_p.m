function Vp = value_p(params,A,F)
        Vp = ((1 - params.beta^params.Tlag) * params.alpha + params.beta^params.Tlag) .* A / (1 - params.beta) - F;
end