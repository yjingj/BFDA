function [hgt, vel, acc] = jolifn(age, coef)
%JOLIFN computes values of Jolicoeur growth model 
%  determined by coefficients in COEF at value in AGE

%  last modified 14 November 2000

a = coef(1);
e = coef(8);
b = coef(2:4);
c = coef(5:7);

n = length(age);
agepe = age + e;
Fvec   = ones(n,1);
DFvec  = zeros(n,1);
D2Fvec = zeros(n,1);
for j=1:3 
   argvec = b(j).*agepe;
   Fvec   = Fvec   +                          argvec.^c(j);
   DFvec  = DFvec  + b(j)   .*c(j).*          argvec.^(c(j)-1);
   D2Fvec = D2Fvec + b(j).^2.*c(j).*(c(j)-1).*argvec.^(c(j)-2);
end
hgt = a.*(1 - 1./Fvec);
vel = (a./Fvec.^2).*DFvec;
acc = -(2.*a./Fvec.^3).*DFvec.^2 + (a./Fvec.^2).*D2Fvec;
