function B = b(params,A,S,k)
   B = params.pI .* (params.rho - 1) .* (A(k) .* sigma(params,S) ./ (params.rho .* params.pI)).^(params.rho ./ (params.rho - 1));
end

