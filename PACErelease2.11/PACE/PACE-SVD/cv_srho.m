function rho_opt = cv_srho(res, x, t_x, y, t_y, nsvd, sig1, method, shrink, regular, rho, ni, tjID)

  mu_x=getVal(getVal(res,'xx'),'mucopy');
  sc_x=getVal(res,'sc_x');
  out1x = getVal(res,'out_x');
  mu_y=getVal(getVal(res,'yy'),'mucopy');
  sc_y=getVal(res,'sc_y');
  out1y = getVal(res,'out_y');
    
  subID = find(ni(:,1) > 1 & ni(:,2) > 1);  %find indices for subjects with ni >= 2

  [muxSub, scxSub] = convertMuPhi(t_x, out1x, mu_x, sc_x, regular);
  [muySub, scySub] = convertMuPhi(t_y, out1y, mu_y, sc_y, regular); 

  cv = zeros(1,length(rho));
  for k = 1:length(rho)
    cv(k) = getSScores2(res, x, t_x,  y, t_y, nsvd, muxSub, scxSub, muySub, scySub, sig1, method, shrink, regular, rho(k,:), subID, tjID);
  end
  IDopt = find(cv == min(cv), 1, 'first');
  rho_opt = rho(IDopt,:);
  fprintf(1, 'Best rho from CV method: (%f, %f) with cv = %f\n', rho_opt(1), rho_opt(2), cv(IDopt));

end