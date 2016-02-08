function R = COR(COV)
% function to generate correlation matrix from the covariance matrix

R = COV ./ (sqrt(diag(COV)) * sqrt(diag(COV))' );

end