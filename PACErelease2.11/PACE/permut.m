function x=permut(a,m)
n=length(a);
for k=n:-1:n-m+1
    temp1=floor(rand(1)*k)+1;
    temp2=a(temp1);
    a(temp1)=a(k);
    a(k)=temp2;
end
x=a(n-m+1:n);
