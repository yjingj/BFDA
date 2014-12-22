% function for computing the GCV bandwidth 
% Input yy: 1 * N vector of observed repeated measurements for all subjects
% Input tt: 1 * N vector of observed time points corresponding to yy
% Input kernel: a character string for the kernel to be used
%               'epan' : epanechikov kernel
%               'gauss' : gaussian kernel
%               'gausvar' : variant of gaussian kernel
%               'rect' : rectangular kernel
%               'quar' : quartic kernel
% Input nwe:    degree of opt polynomial
% Input npoly:  degree of polynomial in the local weighted least square fit
%               ex: linear fit, npoly = 1
% Input nder:   order of derivative in the least square fit, nder =
% 0,1,2,etc 
%               ex: curve fitting nder = 0, first derivative, nder = 1 etc.
% Input regular: if regular = 2, h0 = dstar
%                else,           h0 = 2.5*dstar
% Input verbose: 'on' display diagnostic messages
%                'off' suppress diagnostic messages              
% Input bwmuLocal: 1 local bandwidth choice for mean function
%                  0 global bandwidth choice for mean function
% Output bopt: optimal bandwidth choice among 10 candidates
% Algorithm: 
% 1) find dstar, the minimum distant for every 2+npoly points, ex, dstar = 2+1 = 3 (in linear fit)
%    which guarantees at least three points in a local window. (only for regular = 0)
% 2) Define intial bandwidth h0, which is different for regular = 0,1,2
%    bandwidth candidates are then defined as [h0, h0*q, h0*q^2, ..., h0*q^9]
%    where q = (range(pooled time points)/(4*h0))^(1/9)
% 3) r = range(tt)
%    gcv(h) = 1/(1-(r*k(0))/(Nh))^2 * sum(mu(l)-yy(l))^2
%    here, k(0), kernel evaluated at 0, N : number of pooled observed time points
%          h: one of the candidate bandwidth
%         mu(l): lwls fitted value at time point l, yy(l): observed value
%         at time point l.
% 4) find the h that corresponds to the min(gcv(h))  

function [bopt,gcv] = gcv1d(yy,tt,kernel,nwe,npoly,nder,regular,verbose,bwmuLocal,win)
% min bandwidth to contain at least 2+npoly points in a local window  

r = range(tt);
N = length(tt);  %number of all observed time points
if regular == 0
    dstar = minbwd(tt,2+npoly);
    if dstar > r/4
        dstar = dstar*.75;
        % display(['Warning: the min bandwidth choice is too big, reduce to ' num2str(dstar) '!']);
    end
    h0 = 2.5*dstar;
elseif regular == 1
    h0 = minbwd(tt,1+npoly)*2;
else
    h0 = minbwd(tt,1+npoly)*1.5;
end

if isempty(h0)
    if strcmp(kernel,'gauss') == 1
        h0 = 0.2*r;
    else
        error('The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n');
        bopt = []; gcv = [];
        return;
    end
end

h0 = min(h0, r);

q = (r/(4*h0))^(1/9);
[out1,ignore, idx] = unique(tt);   %distinct pooled observed time points
% fprintf(1,'Bandwidth choices for mean function: \n');
bw = sort(q.^(0:9).*h0) ;                %create 10 h candidates
win1 = ones(1,N);
k0 = kernels(0, kernel);

leave = 0;
nleave = 0;
tooSparse = 0;
while leave == 0
    gcv = Inf*ones(length(bw),1);
    for k = 1:length(bw)
        if length(out1) > 101
            out21 = linspace(min(out1),max(out1),101);
            [invalid, mu] = locpoly(bw(k),kernel, nwe, npoly, nder, tt, yy', win, out21, bwmuLocal);
            if invalid == 0
                newmu = interp1(out21,mu,tt,'spline');
            end
        else
            [invalid, mu] = locpoly(bw(k), kernel, nwe, npoly, nder, tt, yy', win, out1, bwmuLocal);
            if invalid == 0
                newmu = mu(idx);
            end
        end
        % newmu = interp1(out1, mu, tt);
        if invalid == 0
            cvsum = (yy-newmu)*diag(win)*(yy-newmu)'/sum(win);
            gcv(k) = cvsum/(1-(r*k0)/(N*bw(k)))^2;
            if k > 1 && (gcv(k) > gcv(k-1))
                leave = 1;
                break;
            end
        end
    end

    if all(gcv == Inf)
        if nleave == 0 && bw(10) < r
            bopt = bw(10);
            tooSparse = 1;
        else
            error('The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n');
            bopt = [];gcv = [];
            return;
        end
    else
        bopt = bw(find(gcv == min(gcv),1,'first'));
    end

    if bopt == r
        leave = 1;
        % fprintf(1,'Warning: data is too sparse, optimal bandwidth includes all the data! You may want to change to Gaussian kernel!\n');
    elseif bopt == bw(10) && nleave == 0
        if (tooSparse == 1) || (sum(gcv == Inf) == 9)
            % fprintf(1, 'Warning: data is too sparse, retry with larger bandwidths!\n');
            h0 = bw(10)*1.01;
        else
            % fprintf(1, 'Bandwidth candidates are too small, retry with larger choices now!\n');
            h0 = bw(9);
        end
        newr = (0.5:0.05:1).*r;
        id = find(h0 < newr, 1,'first');
        q = (newr(id)/h0)^(1/9);
        bw = q.^(0:9).*h0;              % create 10 h candidates

        bw =  sort(bw) ;
%         if strcmp(verbose, 'on') == 1
%             fprintf(1, 'New bwmu candidates:\n');
%             bw
%         end
    elseif bopt < bw(10) || nleave > 0
        leave = 1;
    end
    nleave = nleave+1;
end 

% if strcmp(kernel,'gauss') == 0 && strcmp(verbose, 'on') == 1
%     fprintf(1,['GCV choice for mean function (npoly = ' num2str(npoly) '): ' num2str(bopt) '\n']);
% end

end
