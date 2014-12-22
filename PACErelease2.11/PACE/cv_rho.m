function rho_opt = cv_rho(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, rho, ni, tjID,verbose)

  subID = find(ni > 1);  %find indices for subjects with ni >= 2

  [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular); 
  LAMBDA = diag(lambda);
  cv = Inf*ones(1,length(rho));

  for k = 1:length(rho)
    cv(k) = getScores2(y, t, mu, phi, lambda, sigma, sig1, noeig, error, method, shrink, out1, regular, muSub, phiSub, LAMBDA, rho(k), subID, tjID);
  end
  IDopt = find(cv == min(cv), 1, 'first');
  rho_opt = rho(IDopt);
  if strcmp(verbose, 'on') == 1
     fprintf(1, 'Best rho from CV method: %f with cv = %f\n', rho_opt, cv(IDopt));
  end
end
