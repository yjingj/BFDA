%This function is mapping the fitted mean function (mu) and
%eigenfunctions (phi) from out1 to the observed times (t)
%[muSub, phiSub] = convertMuPhi(t,out1,mu,phi,regular)
%Input
%    t:  1*n cell array of observed times for n subjects
% out1:  1*N vector of time points for mu and phi
%   mu:  1*N vector of mean vector
%  phi:  N*k matrix of eigenfunctions where k denotes the number
%       of PC
% regular : 0 sparse and irregular data
%           2 completely balance data
%Output:
%   muSub: i) regular = 0
%          1*n cell array of mean functions evaluated at the 
%          observed time points, where muSub{i} is a row
%          vector of mean function evaluated at t{i}
%          ii) regular = 2
%          1*m vector of mean function
%  phiSub: i) regular = 0
%          1*n cell array of eigenfunctions evaluated at the
%          observed time points, where phiSub{i} is a 
%          ni * k matrix of eigenfunctions for ith subject
%          evaluated at t{i}.
%          ii) regular = 2
%          m * k eigenfunctions
function [muSub, phiSub] = convertMuPhi(t, out1, mu, phi, regular)

     if regular == 2
        muSub = mapX1d(out1,mu,t{1});
        phiSub = mapX1d(out1,phi,t{1});
     else
        muSub = cell(1,length(t));
        phiSub = cell(1,length(t));
        for i = 1:length(t)
          muSub{i} = mapX1d(out1,mu,t{i});
          phiSub{i} = mapX1d(out1,phi,t{i});
        end
     end
end
