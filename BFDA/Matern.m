function ac = Matern(d, rho, nu, phi) 
    %d: distance
    %rho: spacial sacle parameter 
    %nu: order of smoothness
    %phi: marginal variance
    if (sum(d < 0) > 0), error('distance argument must be nonnegative'),  end

        dm = d  .* ( sqrt(2 * nu) / rho);
        dm(dm == 0) = 1e-10;
        con = 1 / ((2^(nu - 1)) * gamma(nu));
        ac = phi .* con .* (dm.^nu) .* besselk(nu, dm);
end