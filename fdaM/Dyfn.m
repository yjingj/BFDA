function Dyval = Dyfn(xval, yval, Wfd)
% global Wfd
disp(xval)
Wval  = eval_fd(xval, Wfd);
Dyval = exp(Wval);
