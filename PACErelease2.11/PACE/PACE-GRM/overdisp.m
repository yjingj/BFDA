% Use leave-one-point-out prediction error to choose the overdispersion parameter
% sigma when parametric variance function is used.

function [opt_sigma, opt_pe, odisp, pe]=overdisp(y,t,mux,muy,family,param,lambdax,phix,out1)

    rang = [0.01 10]; 
    odisp = linspace(rang(1),rang(2),50); 
    ncohort = length(y); 
    all = 1:ncohort; 
    pe = [];

    for k = 1:length(odisp)
        [ypred,indpred,yprederr]=grm_predict_lopo(y,t,mux,muy,lambdax,phix,odisp(k),...
             [],family,param,all,out1);
        pe(k)=mean(yprederr);
    end

    opt_pe = min(pe);
    opt_sigma = odisp(pe==opt_pe);
