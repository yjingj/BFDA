% ===========
% Description:
% ===========
% 
%             This program is intended for generalized functional linear regression
%             (GFLM) with known or unknown link and variance functions. The
%             response variable must be a scalar and the predictor is a
%             random trajectory. The predictor trajectories may be measured
%             at irregular times and measurements may be contaminated with
%             noise, in which case PACE will be applied as a preliminary
%             step. Users have the option to specify both the link function and the variance
%             function from a variety of choices provided in the
%             algorithm or to leave either of these functions unspecified. The algorithm
%             uses the Iteratively Weighted Least Squares Method (IWLS) and 
%             SemiParametric Quasi-likelihood Regression (SPQR)
%             See Chiou and M\"uller (2004), M\"uller and Stadtm\"uller
%             (2005) and M\"uller (2005) for details.
% 
%             References: Chiou, J.M., M\"{u}ller, H.G. (2004).  Quasi-likelihood
%             regression with multiple indices and smooth link and variance
%             functions. Scandinavian J. Statistics 31, 367--386.
%             M\"{u}ller, H.G., Stadtm\"uller, U. (2005). Generalized
%             functional linear models. Annals of Statistics 33, 774--805.
%             M\"{u}ller, H.G. (2005). Functional modelling and
%             classification of longitudinal data. Scandinavian J. 
%             Statistics 32, 223--240.
% 
%   
% ======
% Usage:
% ======
% 
%  function [eta,mu,der,sigma2,beta,alpha,tgrid,vbeta,dis,optbw,ggrid,gfct,vgrid,vfct] = iterGFLM(Y,X,T,kernel,bw,datatyp,linktyp,vartyp,theta,fig,pace)
% 
% ======
% Input: 
% ======
%       Y:          1*n vector where Y(i) is the response variable for the
%                   ith subject.
%       X:          n*m matrix or 1*n cell array where X(i,:) 
%                   or X{i} is the predictor for the ith subject.
%       T:          1*m vector or 1*n cell array giving the time
%                   points; if predictors are regular functional
%                   trajectories (X is a matrix), a 1*m vector of the time
%                   points for all subjects; if predictors are irregular
%                   functional trajectories (X is a cell array), a 1*n cell
%                   array where T{i} is the time points for the ith subject.
%       kernel:     kernel function to be used for smoothing link and
%                   variance functions
%                   'epan' - Epanechnikov: 0.75*(1-x^2) on [-1,1] 
%                   'rect' - rectangle:    0.5 on [-1,1]
%                   'gauss'- Gaussian:     exp(-x^2/2)/sqrt(2*pi) on [-4,4]
%                   Default is 'epan'.
%       bw:         a vector of length 3, if specified, will give the bandwidth for smoothing the link function,
%                   the derivative of link function and the variance
%                   function respectively. 
%                   GCV will be used if left unspecified.
%       datatyp:    an integer indicating the type of the predictor.
%                   0 - sparse (or irregular) functional data    
%                   1 - regular data with missing values
%                   2 - completely balanced (regular) data
%                   Default is 0.
%       linktyp:    an integer indicating the type of link function.
%                   0 - unknown; 1 - identity; 2 - power 2; 3 - log; 4 -
%                   logit; 5 - cloglog; 6 - inverse; 7 - sqrt.
%                   Default is 0.
%       vartyp:     an integer indicating the variance function.
%                   0 - unknown; 1 - constant; 2 - binomial; 3 - poisson; 4
%                   - gamma.
%                   Default is 0.
%       theta:      a vector of length n needed only if vartyp is 2 or 4;
%                   theta(i) gives the sample size ni for the ith subject
%                   from binomial B(ni,pi) or gamma G(ni,lambdai)
%                   distribution.
%       fig:        an integer indicating whether the link and variance
%                   functions should be plotted; 0 - no, 1 - yes
%                   Default is 1.
%       pace:       a structure array obtained from setOptions.m that sets parameters for
%                   PACE. Use [] for default choices.
% 
% =======
% Output:  
% =======  
%       eta:       1*n vector where eta(i) is the estimated f(E(Y|mu)) for
%                  the ith subject (here f is the link function). 
%       mu:        1*n vector where mu(i) is the estimated g(eta) (g=f^{-1})
%                  for the ith subject.
%       der:       1*n vector where der(i) is the estimated g'(eta) for the
%                  ith subject.
%       sigma2:    1*n vector where sigma2(i) is the estimated var(Y|mu) for
%                  the ith subject.
%       beta:      a vector giving the regression parameter function,
%                  i.e. eta(i)=alpha+\int beta*X{i} where alpha is the
%                  estimated intercept; if predictors are regular functions
%                  (T is a vector), beta is a 1*m vector evaluated at time
%                  grid T, if predictors are irregularly measured functional
%                  trajectories (T is a 1*n cell array), beta is a 1*100
%                  vector giving the regression function evaluated at time
%                  grid tgrid where tgrid is explained below.
%       alpha:     the estimated intercept (0 if link function is unknown).
%       tgrid:     a 1*m or 1*100 vector giving the evaluated time points
%                  for the estimated regression parameter function beta; if
%                  predictors are regular functional trajectories, tgrid is the same
%                  as T; if predictors are irregular functional
%                  trajectories, tgrid is a 1*100 vector equally spaced in the range
%                  of pooled time points of all trajectories.  
%       vbeta:     covariance matrix of beta; if link function is known, this
%                  is the covariance matrix of [alpha,beta].
%       dis:       the estimated Pearson overdispersion parameter (1 if link
%                  function is unknown).
%       optbw:     1*3 vector giving the final bandwidths for smoothing the
%                  link function, the derivative of link function and the
%                  variance function from generalized cross-validation (GCV); 
%                  NaN if the corresponding function is known. 
%       ggrid:     1*100 vector giving a dense time grid for gfct.
%       gfct:      1*100 vector giving the evaluated link function on time
%                  grid ggrid.
%       vgrid:     1*100 vector giving a dense time grid for vfct.
%       vfct:      1*100 vector giving the evaluated variance function on time
%                  grid vgrid. 

function [eta,mu,der,sigma2,beta,alpha,tgrid,vbeta,dis,optbw,ggrid,gfct,vgrid,vfct] = iterGFLM(Y,X,T,kernel,bw,datatyp,linktyp,vartyp,theta,fig,pace)

if nargin<11 pace=[]; end
if nargin<10|isempty(fig) fig=1; end
if nargin<9 theta=[]; end
if nargin<8|isempty(vartyp) vartyp=0; end
if nargin<7|isempty(linktyp) linktyp=0; end
if nargin<6|isempty(datatyp) datatyp=0; end
if nargin<5 bw=[]; end
if nargin<4|isempty(kernel) kernel='epan'; end

if (linktyp==3|linktyp==7)&length(find(Y<0))>0 error('Some responses are negative!'); end
if (linktyp==4|linktyp==5)&length(find((Y.*(1-Y))<0))>0 error('Some responses are negative or larger than 1!'); end
if linktyp==6&length(find(Y==0))>0 error('Some responses are 0!'); end
if (vartyp==2|vartyp==4)&isempty(theta) error('Theta needs to be given!'); end
if vartyp==2&length(find((Y.*(1-Y))<0))>0 error('Some responses are negative or larger than 1!'); end
if (vartyp==3|vartyp==4)&length(find(Y<0))>0 error('Some responses are negative!'); end

n=length(Y);
% preliminary step: PACE
if iscell(X)
    ngrid=100;
else
    ngrid=size(X,2);
end
if iscell(X)
    tmptmin=min(cell2mat(T));
    tmptmax=max(cell2mat(T));
    newdata=tmptmin:(tmptmax-tmptmin)/(ngrid-1):tmptmax;
    clear tmptmin tmptmax;
    Xnew=X;
    Tnew=T;
else
    newdata=T;
    Xnew=mat2cell(X,ones(1,n),ngrid)';
    Tnew=mat2cell(repmat(T,[n,1]),ones(1,n),ngrid)';
end
if length(pace)==0 pace=setOptions('regular',datatyp,'newdata',newdata,'verbose','off'); end
% fprintf(1,'PACE is running in the preliminary step! \n')
Xpace=FPCA(Xnew,Tnew,pace);
Phi=getVal(Xpace,'phi')';
Xreg=getVal(Xpace,'xi_est')+ones(n,1)*trapz(newdata',Phi'.*(getVal(Xpace,'mu')'*ones(1,size(Phi,1))));
% fprintf(1,['PACE has been applied to the data and the first ' num2str(size(Xreg,2)) ' eigenfunction(s) will be used! \n'])

% initialization
p=size(Xreg,2);
if vartyp==0 | vartyp==1
    sigma2=repmat(var(Y),[1,n]);
else 
    Yadj=iniadj(Y,[],vartyp);
    sigma2=varfun(Yadj,Y,vartyp,theta);
end
if length(bw)==3
    bw12=bw(1:2);
    bw3=bw(3);
else
    bw12=[];
    bw3=[];
end
if linktyp>0
    Yadj=iniadj(Y,linktyp,[]);
    eta=linkfun(Yadj,linktyp,0);
    der=linkfun(Yadj,linktyp,1);
    Xreg2=([ones(1,n); Xreg'])';
    b=betaupd(Xreg2,Y,eta,Y,der,sigma2);
    eta=b*Xreg2';
    mu=linkfun(eta,linktyp,-1);
    der=linkfun(mu,linktyp,1);
    optbw1=NaN;
    optbw2=NaN;
else
    b=repmat(1/sqrt(p),[1,p]); 
    eta=b*Xreg';
    [mu,der,optbw1,optbw2]=gupd(Y,eta,sigma2,kernel,bw12);
end
[sigma2,optbw3]=sigma2upd(Y,theta,mu,vartyp,kernel,bw3);

% iteration
dbeta=1;
iter=0;
% fprintf(1,'Begin iterated reweighted least squares algorithm! \n')
while dbeta>0.01 & iter<50 
    iter=iter+1;
    b_old=b;
    if linktyp>0
        b=betaupd(Xreg2,Y,eta,mu,der,sigma2); 
        eta=b*Xreg2';
        mu=linkfun(eta,linktyp,-1);
        der=linkfun(mu,linktyp,1);
    else
        b=betaupd(Xreg,Y,eta,mu,der,sigma2);
        bnorm=norm(b);
        b=b/bnorm;
        eta=b*Xreg';
        [mu,der,optbw1,optbw2]=gupd(Y,eta,sigma2,kernel,bw12); 
    end 
    [sigma2,optbw3]=sigma2upd(Y,theta,mu,vartyp,kernel,bw3);
    dbeta=norm(b-b_old)/norm(b_old);
end  
% if iter<50
%     fprintf(1,['Solution found after ' num2str(iter) ' iterations ! \n'])
% else
%     fprintf(1,'Maximun number of iterations (50) have been completed. Solution may not converge ! \n')
% end

if linktyp>0
    beta=b(2:end)*Phi;
    alpha=b(1);
else
    beta=b*Phi;
    alpha=0;
end

tgrid=newdata;

if vartyp>0
    if linktyp>0
        dis=sum((Y-mu).^2./sigma2)/(n-p-1);
    else
        dis=sum((Y-mu).^2./sigma2)/(n-p);
    end
else
    dis=1;
end

if vartyp>0
    Vinv=diag((dis*sigma2).^-1);
else
    Vinv=diag(sigma2.^-1);
end
if linktyp>0
    Dx=repmat(der',[1,p+1]).*Xreg2;
    Df=eye(p+1);  
else
    Dx=repmat(der',[1,p]).*Xreg;
    Df=(eye(p)-b'*b)/bnorm;
end
Sigmainv=n*pinv(Dx'*Vinv*Dx);
if linktyp>0
    Phi2=zeros(p+1,size(Phi,2)+1); 
    Phi2(1,1)=1;
    Phi2(2:end,2:end)=Phi;
    vbeta=Phi2'*Df*Sigmainv*Df'*Phi2/n;
else
    vbeta=Phi'*Df*Sigmainv*Df'*Phi/n;
end

optbw=[optbw1,optbw2,optbw3];

ggrid=min(eta):range(eta)/99:max(eta);
vgrid=min(mu):range(mu)/99:max(mu);
if linktyp>0
    gfct=linkfun(ggrid,linktyp,-1);
else
    [invalid,gfct]=locpoly(optbw(1),kernel,[],1,0,eta,mu',1./sigma2,ggrid);
    if invalid==1
        tmptbw=minbwd([eta,ggrid],3);
        clear gfct invalid;
        [invalid,gfct]=locpoly(tmptbw,kernel,[],1,0,eta,mu',1./sigma2,ggrid);
    end
end
if vartyp>1
    vfct=varfun(vgrid,[],vartyp,ones(1,100));
elseif vartyp==1
    vfct=repmat(var(Y-mu),[1,100]);
else
    [invalid2,vfct]=locpoly(optbw(3),kernel,[],1,0,mu,sigma2',ones(1,n),vgrid);
    if invalid2==1
        tmptbw=minbwd([mu,vgrid],3);
        clear vfct invalid2;
        [invalid2,vfct]=locpoly(tmptbw,kernel,[],1,0,mu,sigma2',ones(1,n),vgrid);
    end
end
if fig==1
    figure
    plot(eta,Y,'.r')
    hold on
    plot(ggrid,gfct,'-b')
    hold off
    xlabel('estimated eta')
    ylabel('observed Y')
    if linktyp==0 
        title('estimated link function') 
    else
        title('link function') 
    end
    figure
    plot(mu,(Y-mu).^2,'.r')
    hold on
    plot(vgrid,dis*vfct,'-b')
    hold off
    xlabel('estimated mu')
    ylabel('estimated squared residual')
    if vartyp==0 
        title('estimated variance function') 
    else
        title('variance function with overdispersion') 
    end
end

end
    