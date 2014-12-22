%function [raw_ccov, tpairn] = getRawCCOV(x,t_x, y,t_y, mu_x, mu_y, isYFun, regular)
%obtain the raw cross covariance surface.
%Input:
% x    :  1*n cell array for predictor function x,  x{i} is the vector of measurements for the ith subject,
%         i=1,...,n.
%         It contains subjects that are used for prediction, when "isNewSub" is defined to be case a) or b),
%         that is, nn number of subjects will be used for estimation and the n-nn subjects will be used for
%         prediction.
% t_x  :  1*n cell array, t_x{i} is the vector of time points for the ith subject on which
%         corresponding measurements x{i} are taken, i=1,...,n.
%         It contains subjects that are used for prediction, when "isNewSub" is defined to be case a) or b),
%         that is,nn number of subjects will be used for estimation and the n-nn subjects will be used for
%         prediction.
% y    :  1*n cell array for response function y,  y{i} is the response vector for the ith subject,
%         i = 1,..., n
%         or 1*n vector for scalar response y, y(i) is the response value for the ith subject,
%         i = 1,...,n
% t_y  :  1*n cell array, t_y{i} is the vector of time points for the ith subject on which
%         corresponding measurements y{i} are taken, i = 1,...,n.
%         or [] if y is scalar response.
% mu_x :  1*N vector, estimated mean functions valued at distinct input time ponits with 
%         ascending order from all subjects for x.
% mu_y :  1*N vector, estimated mean functions valued at distinct input time ponits with 
%         ascending order from all subjects for y.
% isYFun: a logical value, isYFun = 1 (Y is a functional response)
%                          isYFun = 0 (Y is a scalar response)
% regular:0, sparse functional data.
%         1, irregular but not sparse grid data, allowing some missing or different time points.
%         2, completely balanced data.
%Output:
% tpairn  : For functional Y (isYFun = 1), 2*NN vector denotes the pairs of time points for subject 
%           concatenating as two vectors
%           For scalar Y (isYFun = 0), 1*NN vector denotes the time points for X
% raw_ccov: 1 * NN vector of raw covariance corresponding to tpairn 

function [raw_ccov, tpairn] = getRawCCOV(x,t_x, y,t_y, mu_x, mu_y, isYFun, regular)
tt_x = []; tt_y = [];
xx = []; yy = [];
n = length(x);

if isYFun==1                                   %Y is functional
    if regular ~= 2                              %sparse irregular or regular with missing data
        for i = 1:length(t_x)
            [tx,ty] = meshgrid(t_x{i},t_y{i});
            [x1,y1] = meshgrid(x{i},y{i});
            tt_x = [tt_x; tx(:)];                    %vectorize observed time points
            tt_y = [tt_y; ty(:)];                    %vectorize observed time points
            xx = [xx; x1(:)];                        %vectorize observed measurements
            yy = [yy; y1(:)];                        %vectorize observed measurements
        end
        [outgrid1, ignore, id1] = unique(tt_x);
        [outgrid2, ignore, id2] = unique(tt_y);
        raw_ccov = (xx'-mu_x(id1)).*(yy'-mu_y(id2));
        tpairn = [tt_x tt_y]';
    else                                             %regular balance data
        xx = reshape(cell2mat(x),length(x{1}), n)';
        mux = repmat(mu_x,n,1);
        yy = reshape(cell2mat(y),length(y{1}), n)';
        muy = repmat(mu_y,n,1);
        raw_ccov = (xx-mux)'*(yy-muy)/n;
        raw_ccov = raw_ccov(:)';
        [outtx,outty] = meshgrid(t_x{1},t_y{1});
        tpairn = [outtx(:)'; outty(:)'];
    end
else %When y is a scalar, the raw covariance is a function of T_{ij} for X
    if regular ~= 2                              %sparse irregular or regular with missing data
        for i = 1:length(t_x)
            tt_x = [tt_x; t_x{i}(:)];                    %vectorize observed time points
            xx = [xx; x{i}(:)];                          %vectorize observed measurements
            yy = [yy; repmat(y(i),length(t_x{i}),1)];    %vectorize observed measurements
        end
        [outgrid1, ignore, id1] = unique(tt_x);
        raw_ccov = (xx'-mu_x(id1)).*(yy'-repmat(mu_y,1,length(yy)));
        tpairn = tt_x';
    else                                             %regular balance data
        xx = reshape(cell2mat(x),length(x{1}), n)';
        mux = repmat(mu_x,n,1);
        muy = repmat(mu_y,1,n);
        raw_ccov = (y-muy)*(xx-mux)/n;
        tpairn = t_x{1};
    end
end
end
