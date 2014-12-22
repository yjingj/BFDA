function [ypred, rand_int, ypred_int]=pred_int(y,t,mu,phi,out1,omit,ind_pred,xi_est,xi_var,tout);
ypred={};
rand_int={};
ypred_int={};
for i=1:length(ind_pred)
    get=ind_pred(i);
    phii=[];
    mutemp=interp1(out1,mu,tout{get});
    for k=1:size(phi,2)
        phii(:,k)=interp1(out1,phi(:,k),tout{get})';
    end
    ypred{get}=mutemp+xi_est(i,:)*phii';
    rand_int = [];
    ypred_int = [];
 %   var0=diag(phii*xi_var{get}*phii');
 %   rand_int{get}=[xi_est(i,:)*phii'-1.96*sqrt(var0');...
 %                  xi_est(i,:)*phii'+1.96*sqrt(var0')];
 %   ypred_int{get}=[mutemp+rand_int{get}(1,:);mutemp+rand_int{get}(2,:)];
end
