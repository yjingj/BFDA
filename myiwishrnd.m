function [iK, K]  = myiwishrnd(df, G, p)

  iG = pinv(G);
  L = mychol(iG);
  %try L = chol(iG, 'lower'); catch L =  chol(topdm(iG), 'lower'); end
  X = L * normrnd(0, 1, p, df);
  K = X * X';
  iK = pinv(K); 

end