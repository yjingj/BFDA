% This function is to estimate the mean functions of the observed non-Gaussian
% process Y and the latent Gaussian process X, evaluated at output grid
% out1.
function [muy,bw_muy,mux,out1]=grm_mu_NP(y,t,family,param,kernel,bwmuy,bwmuy_gcv,regular,verbose)

    tt = cell2mat(t);  % 1 x N vector to hold the observed time points from all subjects
    yy = cell2mat(y);  % 1 x N vector to hold the observed measurements from all subjects
    out1=unique(sort(tt));

    if bwmuy == 0
        if bwmuy_gcv == 1
            %use GCV method to choose bw for mean function
            bw_muy = gcv_lwls(yy,tt,kernel,1,1,0,regular,verbose); 
            bw_muy = adjustBW1(kernel,bw_muy,1,0,regular,verbose);
        else
            %use CV method to choose bw for mean function
            bw_muy=cvfda_lwls(y,t,kernel,1,1,0,regular,verbose);
        end
    elseif bwmuy > 0
        bw_muy = bwmuy;
    else
        fprintf(1,'Error: Bandwidth choice for the mean function must be positive!\n');
        return;
    end

    %define the vector of case weight in the local weighted least square
    %here, it is set to be one for all subjects
    win1=ones(1,length(tt));
    [invalid, muy] = lwls(bw_muy,kernel,1,1,0,tt,yy',win1,out1);

    mux = glink_inv(muy,family,param);
