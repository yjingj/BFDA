%function [bw_ccov,bw_cddt, gcv]=gcv2_mullwlsn(t,ngrid,regular,kernel,rcov)
%This function computes the bandwidth choices for the smooth covariance surface
%using GCV method by pooling all time point pairs from all subjects together.
function [bw_ccov, bw_cddt, gcv]=gcv2_mullwlsn(t_x,t_y,ngrid,regular,kernel,tpair,cxx,win)
ttx = cell2mat(t_x);
tty = cell2mat(t_y);
out1x=unique(ttx);
out1y=unique(tty);
clear ttx tty;
a0x=min(out1x);
b0x=max(out1x);
a0y=min(out1y);
b0y=max(out1y);

h0x = getMinb(t_x,out1x,regular);
h0y = getMinb(t_y,out1y,regular);
if strcmp(kernel,'gauss')
    h0x = h0x*0.2;
    h0y = h0y*0.2;
end
%t1 = cputime;
% Calculate the raw covariance for various cases
% depending on regular, error
%rcov = getRawCov(y,t,mu,regular, error);
clear t_x t_y;

%fprintf(1,['Time after the raw covariance calculation: ' num2str(cputime-t1) '\n']);
% create 10 bandwidth candidates

%Alternative choice: h0= 2.5*dstar; or h0 = 1.1*dstar;
%h0 = dstar;
N = length(cxx);   %No. of pairs for the covariance
rx = range(out1x);
ry = range(out1y);
clear out1x out1y;

%Alternative way to set h0
%dstar = minb(tpairn(1,:),10);
%h0 = sqrt(2.5*dstar/range(tpairn(1,:)));
qx = (rx/(4*h0x))^(1/9);
bwx = qx.^(0:9).*h0x;        %create 10 h candidates
qy = (ry/(4*h0y))^(1/9);
bwy = qy.^(0:9).*h0y;        %create 10 h candidates
%fprintf(1, 'bwccov candidates:\n')

bw = [bwx' bwy'];     %same bandwidth choice for each direction of the covariance

%bw

k0 = mykernel(0, kernel);

out2x = linspace(a0x,b0x,ngrid); %then length(out2x) = ngrid
out2y = linspace(a0y,b0y,ngrid); %then length(out2y) = ngrid
bw = sort(bw);
%curTIME = cputime;

leave = 0;
nleave = 0;
while leave == 0
    gcv = Inf*ones(size(bw,1),1);
    for k = 1:size(bw,1)

        %    fprintf(1,['k = ' num2str(k) '\n']);
        if regular == 1
            [invalid, ccov]= mullwlsk(bw(k,:), kernel, tpair, cxx', win, out2x, out2y,rcovcount);
        else
            [invalid, ccov]= mullwlsk(bw(k,:), kernel, tpair, cxx', win, out2x, out2y);
        end
        %      fprintf(1,'Done with smooth the raw covariance\n');
        %interpolate the smooth covariance from (out2x,out2y) to (tpairn(1,:), tpairn(2,:))
        if invalid ~= 1
            newccov = interp2(out2x,out2y,ccov, tpair(1,:),tpair(2,:),'spline');
            clear ccov;
            %      fprintf(1,'Done with the interpolation\n');
            cvsum = (cxx-newccov)*(cxx-newccov)';
            clear newccov;
            bottom = 1-(1/N)*((rx*ry*k0^2)/bw(k,1)/bw(k,2));
            gcv(k) = cvsum/(bottom)^2;
            tmp = gcv(~isinf(gcv));
            if length(tmp) > 1 && gcv(k) > gcv(k-1)
                leave = 1;
                break;
            end
        end
    end
    if all(gcv == Inf)
        fprintf(1,'Warning: the data is too sparse!\n');
        bw_ccov = max(bw)*2;
    else
        bw_ccov = bw(find(gcv == min(gcv),1,'first'),:);
    end

    if bw_ccov(1) == bw(10,1) && nleave == 0
        fprintf(1, 'Bandwidth candidates are too small, retry with larger choices now!\n');
        h0x = bw(9,1);
        h0y = bw(9,2);
        qx = (rx/(4*h0x))^(1/9);
        bwx = qx.^(0:9).*h0x;        %create 10 h candidates
        qy = (ry/(4*h0y))^(1/9);
        bwy = qy.^(0:9).*h0y;        %create 10 h candidates
        fprintf(1, 'New bwccov candidates:\n');
        bw = [bwx' bwy'];     %same bandwidth choice for each direction of the covariance
        bw
    elseif bw_ccov(1) == bw(1,1) && nleave == 0
        fprintf(1, 'Bandwidth candidates are too large, retry with smaller choices now!\n');
        h0x = bw(2,1)/4;
        h0y = bw(2,2)/4;
        qx = (4)^(1/9);
        bwx = qx.^(0:9).*h0x;        %create 10 h candidates
        qy = (4)^(1/9);
        bwy = qy.^(0:9).*h0y;        %create 10 h candidates
        fprintf(1, 'New bwccov candidates:\n');
        bw = [bwx' bwy'];     %same bandwidth choice for each direction of the covariance
        bw
    elseif bw_ccov(1) < bw(10,1) || nleave > 0
        leave = 1;
    end
    nleave = nleave+1;

end  %end of while loop

bw_cddt = bw;
% if strcmp(kernel,'gauss')
%     bw_ccov = bw_ccov*1.1;
% end

fprintf(1,['GCV bandwidth choice for COV function: (' num2str(bw_ccov(1)) ',' num2str(bw_ccov(2)) ')\n']);

end
