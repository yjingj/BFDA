function [L] = mychol(sig) 

[U, S, V] = svd(sig);
L = U * sqrt(S);

end 

