function p = Covprob(x, xl, xu, dim)

switch dim
    case 1
        p = mean((x >= xl) & (x <= xu));
    case 2
        p = mean(mean((x >= xl) & (x <= xu)));
end