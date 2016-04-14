function C = mtimesw(A, B, W)
%  generalizes mtimes or * to allow for a metric vector or matrix
%  If W is a matrix, the function computes A*W*B
%  if W is a vector, the function computes A*diag(W)*B

[arow, acol] = size(A);
[brow, bcol] = size(B);
[wrow, wcol] = size(W);

if acol ~= brow
    error('no. cols of A not equal to no. rows of B');
end

if wcol == 1 
    if wrow ~= acol
        error('length of W not equal to no. rows of B');
    end
    W = repmat(W,1,bcol);
    C = A*(W.*B);
    return;
end

if wrow == 1 
    if wcol ~= acol
        error('length of W not equal to no. rows of B');
    end
    W = repmat(W,bcol,1);
    C = A*(W.*B);
    return;
end

if wcol ~= acol || wrow ~= acol
    error('W is a matrix but of wrong dimensions.');
end

C = A*W*B;


    