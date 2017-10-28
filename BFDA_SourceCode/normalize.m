function [output] = normalize(X)
% normalize columns of X
    [~, n] = size(X);
    for i = 1:n 
        X(:, i) = X(:, i) - mean(X(:, i));
        X(:, i) = X(:, i) - std(X(:, i));
    end
    output = X;
end
