% function for updating mu and der in the iteration
% INPUTs: Y,eta,sigma2,kernel as in the inputs of iterGFLM; bw is a 1*2 vector
% giving the bandwidths for smoothing link function and the derivative of
% link function
% OUTPUTs: mu,der as in the outputs of iterGFLM; optbw1,optbw2 give the GCV
% (if not specified) bandwidths for smoothing link function and the
% derivative of link function

function [mu,der,optbw1,optbw2] = gupd(Y,eta,sigma2,kernel,bw)

% for sigma2=0, use the the nonzero minimum instead
tmptindx=find(sigma2>0);
if length(tmptindx)==0
    win=repmat(1,[1,length(Y)]);
else
    win=repmat(max(1./sigma2(tmptindx)),[1,length(Y)]);
end
win(tmptindx)=1./sigma2(tmptindx);
clear tmptindx
if length(bw)==2
    optbw1=bw(1);
    optbw2=bw(2);
else
    optbw1=gcv1d(Y,eta,kernel,[],1,0,2,'off',0,win);
    optbw2=gcv1d(Y,eta,kernel,[],2,1,2,'off',0,win);
end
[tmpt,mu]=locpoly(optbw1,kernel,[],1,0,eta,Y',win,eta);
clear tmpt;
[tmpt,der]=locpoly(optbw2,kernel,[],2,1,eta,Y',win,eta);
clear tmpt;

end