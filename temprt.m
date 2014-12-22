%%
%%%%BHM  
rho = sqrt(2 * nu) / rt;
    A = Matern(D, rho, nu, 1); % update A  
    
%update gamma_theta1 (1/theta1)  
    lnsd = 0.1;
    rt1 = lognrnd(log(rt), lnsd);
    rtP(iter) = rt1;
    ratio = min(fTheta1(rt1, rs, D, nu, K, at, bt, m)- fTheta1(rt, rs, D, nu, K, at, bt, m)...
        - log(lognpdf(rt1, log(rt), lnsd)) + log(lognpdf(rt, log(rt1), lnsd)), 0);
    if(log(rand(1)) < ratio)
      rt = rt1;  
      count = count +1;
    end
    
    
  %%
  %%%% Matern model
     %rs = gamrnd(p*n/2 + as, ...
     %  1/(bs + trace( (Z - repmat(mu, 1, n))' * pinv(A) * (Z - repmat(mu, 1, n)) ) / 2 ));
     
   %update rt
    lnsd = 0.1;
    rt1 = lognrnd(log(rt), lnsd);
    rtP(iter) = rt1;
    ratio = min(fTheta1_matern(rt1, rs, D, nu, Z, mu, at, bt, n) - ...
               fTheta1_matern(rt, rs, D, nu, Z, mu, at, bt, n)...
        - log(lognpdf(rt1, log(rt), lnsd)) + log(lognpdf(rt, log(rt1), lnsd)) , 0);
    if(log(rand(1)) < ratio)
      rt = rt1;  
      count = count +1;
    end