function mu=mu_true(t,p)

  t(~(t >= 0 & t <= p)) = 0;
  mu = t+sin(t);

end
