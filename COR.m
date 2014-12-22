function R = COR(COV)
   
R = COV ./ (sqrt(diag(COV)) * sqrt(diag(COV))' );

end