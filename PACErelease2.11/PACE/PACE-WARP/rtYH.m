function ytrans=rtYH(t,y,h)
% evaluates functional curve (t,y) at different time points h
% input 
% t, 1*m or m*1 vector, gives the time points
% y, n*m matrix or n*1 struct where each row gives a function evaluated at time t or a spline 
% h, a vector of arbitrary size which gives evaluated time points
% output
% ytrans, n*m matrix where each row gives the evaluated values at time h

n=size(y,1);
ytrans=[];
if 1-isnumeric(y)
    yfd=y;
else
    for i=1:n
        yfd(i,:)=spapi(2,t,y(i,:));
    end
end         
for i=1:n   
    ytrans(i,:)=fnval(yfd(i,:),h(i,:));
end

end