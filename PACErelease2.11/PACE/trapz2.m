function intg=trapz2(z,x,y);
temp=[];
n1=length(x);
for i=1:n1
    temp(i)=trapz(y,z(i,:));
end
intg=trapz(x,temp);
