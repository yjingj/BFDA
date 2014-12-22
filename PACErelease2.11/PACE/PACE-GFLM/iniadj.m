% function to addjust Y for initial estimates of mu, eta or sigma2 

function Yadj = iniadj(Y,linktyp,vartyp)

if (linktyp==3 | linktyp==7) & length(find(Y==0))>0
    Yadj=Y;
    Yadj(find(Y==0))=Yadj(find(Y==0))+min([Y(find(Y>0)),0.01]);
elseif (linktyp==4 | linktyp==5) & length(find((Y.*(1-Y))==0))>0
    tmptindx=find(Y.*(1-Y)==0);
    Yadj=Y;
    Yadj(tmptindx)=Yadj(tmptindx)+sign(0.5-Yadj(tmptindx)).*min([Y(find(Y>0)),1-Y(find(Y<1)),0.01]);
elseif linktyp==6 & length(find(Y==0))>0
    Yadj=Y;
    Yadj(find(Y==0))=Yadj(find(Y==0))+sign(rand(1,lengt(find(Y==0)))-0.5)*hmin([Y(find(abs(Y)>0)),0.01]);
elseif length(linktyp)==1
    Yadj=Y;
end

if vartyp==2 & length(find((Y.*(1-Y))==0))>0
    tmptindx=find(Y.*(1-Y)==0);
    Yadj=Y;
    Yadj(tmptindx)=Yadj(tmptindx)+sign(0.5-Yadj(tmptindx)).*min([Y(find(Y>0)),1-Y(find(Y<1)),0.01]);
elseif length(vartyp)==1 & length(find(Y==0))>0
    Yadj=Y;
    tmptindx=find(Y==0);
    Yadj(tmptindx)=ones(1,length(tmptindx))*min([Y(find(Y>0)),0.01]);
elseif length(vartyp)==1
    Yadj=Y;
end
    
end