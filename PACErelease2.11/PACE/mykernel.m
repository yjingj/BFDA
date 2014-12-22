%[kx] = mykernel(x,kernel)
%This function calculates the kernel density based on the
%input vector x
%Input x: a vector
%Input kernel: a chacracter string for one the following:
%'quar' : quartic kernel
%'epan' : epanechikov kernel
%'gauss' : gaussian kernel
%'gausvar' : variant of gaussian kernel
%'rect' : rectangular kernel
%Output kx: a vector of the kernel density of x
function [kx] =  mykernel(x, kernel)
     if strcmp(kernel,'quar')
        kx = (15/16).*(abs(x) <= 1).*(1-x.*x).^2;
     elseif strcmp(kernel,'epan')
        kx = 0.75.*(abs(x) <= 1).*(1-x.*x);
     elseif strcmp(kernel,'rect')
        kx = 0.5.*(abs(x) <= 1);
     elseif strcmp(kernel,'gausvar')
        kx = (1/sqrt(2*pi)).*exp(-0.5.*x.*x).*(1.25-0.25*x.^2);
     elseif strcmp(kernel, 'gausvar1')
	kx = (1/sqrt(2*pi)).*exp(-0.5.*x.*x).*(1.5-0.5*x.^2);
     elseif strcmp(kernel,'gausvar2')
        k1 = (1/sqrt(2*pi)).*exp(-0.5.*x.*x).*(1.25-0.25*x.^2);
        k2 = (1/sqrt(2*pi)).*exp(-0.5.*x.*x).*(1.5-0.5*x.^2);
        kx = k1.*k2;
     else
        kx = (1/sqrt(2*pi)).*exp(-0.5.*x.*x);
     end
