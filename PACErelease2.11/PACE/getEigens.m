function [lambda, phi, eigen, noeig] = getEigens(xcov,out1,out21,noeig, varargin)

h=range(out21)/(length(out21)-1);
%calculate the first noopt eigenfunctions and eigenvalues
%[u d eigen]=svds(xcov,noeig,'L');  %it always results in positive eigenvalues
opts.disp = 0;  %don't show the intermediate steps
ngrid = size(xcov,1);

[eigen d] = eigs(xcov,ngrid-2,'lm',opts);

d = diag(d);
iddx = find(imag(d));                           %remove any imaginary eigenvalues                                               
if isempty(iddx) == 0 
   %if length(varargin) > 0                                                                                                    
      %fprintf(1,['Warning: ' num2str(length(iddx)) ' eigenvalues are complex and are removed!\n']);
   %end
   error('Error: %d eigenvalues are complex. The estimated auto-covariance surface is not symmetric!',length(iddx));
   %iddx = find(~image(d));
   %d = d(iddx);
   %eigen = eigen(:,iddx);                       
end                                                                                                                      
iddx = find(d <= 0);                                                                                                       
if isempty(iddx) == 0   
   if length(varargin) > 0                                                                                                     
     fprintf(1,['Warning: ' num2str(length(iddx)) ' real eigenvalues are negative or zero and are removed!\n']);
   end                    
   eigen = eigen(:,d > 0);                                  
   d = d(d > 0);                                  %retain only the positive eigenvalues
end                                              

if noeig > length(d)
    noeig = length(d);
    fprintf(1,['Warning: at most ' num2str(noeig) ' number of PC can be selected!\n']);
end

eigen = eigen(:,1:noeig);
d = d(1:noeig);   

eigen = eigen./sqrt(h);
%lambda = h*diag(d)';
lambda = h*d';

%for i=1:noeig
%  eigenNorm(i) = trapz(out21, eigen(:,i).^2);
%  lambda(i)=d(i,i)*eigenNorm(i);
%  eigen(:,i)=v(:,i)./sqrt(h);
%end


%Check for the validity of eigen functions
%remove any non-positive lambda and update the
%corresponding "noeig" -- no. of PC
%negLambda = find(lambda <= 0);
%if isempty(negLambda) == 0
%   for jj = 1:length(negLambda)
%     aa = lambda(negLambda(jj));
%     display(['Warning: lambda(' num2str(negLambda(jj)) ') = ' num2str(aa) ' is either negative or zero and is removed!']);
%   end
%     posLambda = find(lambda > 0);
%     noeig = length(posLambda);
%     fprintf(1,['The resulting number of PC : ' num2str(noeig) '\n']);
%     lambda = lambda(posLambda);
%     eigen = eigen(:, posLambda);
%end



%normalized smoothed eigen functions
%This eigen function is without interpolation and
%evaluated at out21

%ww = out21 - [0 out21(1:length(out21)-1)];
%ttemp = (1:length(ww))'; 
%WW = sparse(ttemp, ttemp, ww');

for i = 1:noeig
    %eigen(:,i) = eigen(:,i)/sqrt(eigen(:,i)'*WW*eigen(:,i));        
    %eigen(:,i) = eigen(:,i)/sqrt(romb2(out21,eigen(:,i).^2));
    eigen(:,i) = eigen(:,i)/sqrt(trapz(out21,eigen(:,i).^2));    
    if eigen(2,i)< eigen(1,i)
       eigen(:,i)=-eigen(:,i);
    end
end

%w=out1-[0,out1(1:length(out1)-1)];
%temp=(1:length(w))';
%W=sparse(temp,temp,w');

% interpolate from the normalized the eigenfunctions
phi=interp1(out21, eigen, out1, 'spline');
if noeig == 1
   phi = phi';
end
% normalize smoothed eigenfunctions
for i=1:noeig
   %phi(:,i) = phi(:,i)/sqrt(phi(:,i)'*W*phi(:,i));
   %phi(:,i) = phi(:,i)/sqrt(romb2(out1,phi(:,i).^2));
   phi(:,i) = phi(:,i)/sqrt(trapz(out1,phi(:,i).^2));
   %if phi(2,i)<phi(1,i)
   %   phi(:,i)=-phi(:,i);
   %end
end
