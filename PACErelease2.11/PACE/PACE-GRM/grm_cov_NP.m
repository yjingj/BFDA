% This function is to estimate the covariance surfaces of the observed non-Gaussian
% process Y and the latent Gaussian process X, evaluated at output grid
% out21.
function [invalid,covy, bw_covy, covx, out21]=grm_cov_NP(y,t,muy,mux,family,param,...
    kernel,bwcovy,bwcovy_gcv,ntest1,ngrid1,ngrid,regular,error,out1,verbose)

rcov = getRawCov(y,t,out1,muy,regular,0);           %obtain raw covariance;

% Choose bandwidth for covariance estimation
if bwcovy(1)==0 || bwcovy(2)==0
    if bwcovy_gcv == 1
        [bw_covy] = gcv_mullwlsn(t,ngrid1,regular,error,kernel,rcov,verbose);
        bw_covy = adjustBW2(kernel,bw_covy,1,0,regular,verbose);
    else
        [bw_covy] = cv_mullwlsn(y,t,muy,ntest1,ngrid1,regular,error,kernel,rcov,verbose);
    end
elseif bwcovy > 0
    bw_covy = bwcovy;
elseif bwcovy(1) < 0 || bwcovy(2) < 0
    fprintf(1,'Error: Bandwidth choice for the covariance function must be positive!\n');
    return;
end

% Estimate the covariance with the selected or defined bandwidth
out21 = linspace(min(out1),max(out1),ngrid);
rcov1 = rcov;
if error == 1
    tpairn = rcov1.tpairn;
    tneq = find(tpairn(1,:)~=tpairn(2,:));
    cyy = rcov1.cyy;
    rcov1.tpairn = tpairn(:,tneq);
    rcov1.cxxn = cyy(tneq);
    rcov1.win = ones(1,length(rcov1.cxxn));
    if regular == 1
        rcov1.count = rcov1.count(tneq);
    end
end

if regular == 1
    [invalid,covy]=mullwlsk(bw_covy,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21,rcov1.count);  %smooth raw covariance;
else
    [invalid,covy]=mullwlsk(bw_covy,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21);  %smooth raw covariance;
end
covy = (covy+covy')/2;

% Obtain the original and nonnegative-definitely adjusted covariance estimates 
% for the hidden process X
if invalid == 0
    mux1 = interp1(out1,mux,out21,'spline');
    covx = diag(1./glink_der(mux1,family,param))*covy*diag(1./glink_der(mux1,family,param));
    covx = (covx+covx')/2;
%     [v,d] = eig(covx); 
%     d(d<0) = 0;
%     d = diag(d);
%     covxp = v*diag(d)*v';
elseif invalid == 1
    covx=[]; % covxp=[];
    fprintf(1,'Error: Bandwidth choice for the covariance function should be increased!\n');
    return;
end
