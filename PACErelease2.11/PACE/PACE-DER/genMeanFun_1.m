function [mu] = genMeanFun_1(t,der,sigma)
     if nargin == 1
       der = 0;
     end
     if nargin <= 2
       sigma = 0.1;
     end
     mu = 0.5;
     if der == 0
       %mu = 3-3*t+exp(-((t-0.5).^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);%decreasing
        mu = 4*t+exp(-((t-0.5).^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);  %increasing
     elseif der == 1
       %mu = -3-((t-mu)/(sqrt(2*pi)*sigma^3)).*exp(-((t-0.5).^2)/(2*sigma^2));
        mu = 4-((t-mu)/(sqrt(2*pi)*sigma^3)).*exp(-((t-0.5).^2)/(2*sigma^2));
     elseif der == 2
        mu = (1/(sqrt(2*pi)*sigma^3))*(((t-mu)/sigma).^2-1).*exp(-((t-0.5).^2)/(2*sigma^2));
     end
    
end
