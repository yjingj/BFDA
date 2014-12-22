function qdm=qdmscores1(xi)
% Arrange the linear and crossterms of FPC scores to form a whole 
% design matrix for the functional quadratic model

[n,K]=size(xi); qdm=[];

temp1=xi;
temp2=[];
for k=1:K
    for j=1:k
        temp2=[temp2,xi(:,k).*xi(:,j)];
    end
end
qdm=[temp1,temp2];
        
