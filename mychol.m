function [L] = mychol(sig) 

[U, S, V] = svd(sig);
DS = diag(S);
DS = sqrt(DS);
minDS = mean(DS)*10^(-8);
DS(DS < 0) = minDS;
L = U * diag(DS);

end 

