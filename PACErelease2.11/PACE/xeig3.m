function phi=xeig3(t,T)
phi=zeros(3,length(t));
if t>=0 & t<=T
    phi(1,:)=-sqrt(2/T)*cos(2*pi*t/T);
    phi(2,:)=sqrt(2/T)*sin(2*pi*t/T);
    phi(3,:)=-sqrt(2/T)*cos(4*pi*t/T);
else
    fprintf(1,'Error: t should be smaller than or euqual to T!\n');
    phi = [];
end      
