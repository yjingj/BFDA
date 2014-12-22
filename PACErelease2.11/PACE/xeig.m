%Create eigenfunctions in the domain of [0,T]
%function phi = xeig(t,T,numPhi)
%Input:
%    t: 1*m vector of time points for each of the eigenfunctions
%    T: a real value number where T > 0, the data is to be 
%       sampled in the domain of [0,T] 
%  numPhi: positive integer for the number of eigenfunctions to be returned
%Output:
%   phi: numPhi * m matrix of eigenfunctions
%        phi(k,:) is the kth eigenfunction
% ex:
%   t = 0:0.01:1;
%   T = 1;
%   phi1 = xeig(t,T,1)  %1 eigenfunction
%   phi2 = xeig(t,T,2)  %2 eigenfunctions
%   phi3 = xeig(t,T,3)  %3 eigenfunctions
function phi=xeig(t,T,numPhi)

  if nargin == 2
   numPhi = 2;
  end

  id = find(t < 0 | t > T);
  if isempty(id) == 0
    fprintf(1,['Warning: t must be in [0, ' num2str(T) ']. Invalid t''s are removed!\n']);
    t(id) = [];
  end

  if numPhi == 1
    phi = -sqrt(2/T)*cos(2*pi*t/T);
    return;
  end 

  if numPhi == 2
      phi(1,:) = -sqrt(2/T)*cos(2*pi*t/T);
      phi(2,:) = sqrt(2/T)*sin(2*pi*t/T);
      return;
  end



  phi=zeros(numPhi,length(t));
  id = 1:numPhi;
  oddID = mod(id, 2);
  oddFactor = 1:sum(oddID);
  evenID = oddID == 0;
  evenFactor = 1:sum(evenID);

  phiOdd = zeros(sum(oddID),length(t));
  phiEven = zeros(sum(evenID), length(t));

  for i = 1:sum(oddID)
     phiOdd(i,:) = -sqrt(2/T)*cos(2*oddFactor(i)*pi*t/T);
  end
  phi(find(oddID == 1),:)= phiOdd;

  for i = 1:sum(evenID)
     phiEven(i,:) = sqrt(2/T)*sin(2*evenFactor(i)*pi*t/T);
  end
  phi(find(evenID== 1),:)= phiEven;

end

