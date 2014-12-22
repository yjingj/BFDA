function [phi_sign] = checkPhiSign(out1, phi, true_phi)

  a1 = trapz(out1, (phi-true_phi).^2);
  a2 = trapz(out1, (-phi-true_phi).^2);

  if a1 <= a2
    phi_sign = 1;
  else
    phi_sign = -1;
  end

end
