function a = area(x, t)
    % calculate area of the functional curve x(t)
    
    tdiff = diff(t);
    n = length(x);
    if n>1
    a =  0.5 * sum ( abs(x(2:n) + x(1:(n-1)) ) .* tdiff) ;
    else
        error('length of x should > 1')
    end

end