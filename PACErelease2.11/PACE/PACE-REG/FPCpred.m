%function [newy newcb] = FPCpred(res, xx, yy, newx, new_tx, y, new_ty)
%======
%Input:
%======
%  res, xx, yy:  The returned values from FPCreg. See FPCreg() for more details.
%  newx:   1*numNewSub cell array contains measurements for new x functions
%  new_tx: 1*numNewSub cell array contains time points corresponding to the newx
%  y:      the scalar response Y vector used for estimation 
%          Only needed when Y is a scalar.
%  new_ty: 1*numNewSub cell array contains time points.
%  FIT:    an indicator with values 0 or -1.
%          Refer to the input FIT of FPCreg for more details.
%  alpha:  the level of the confidence bands.  alpha = 0.05 if is left empty. 
%          No confidence bands will be produced if the inputted alpha is out of (0, 1).
%
%=======
%Output:
%=======
%  newy:   1*numNewSub cell array contains predicted measurements for the
%          corresponding newx and new_ty
%  newcb:  when y is functional, 2*numNewSub cell array, newcb{1,i} is the
%          lower confidence bound for the i-th new subject, while newcb{2,i} 
%          is the upper confidence bound;
%          when y is scalar, 2*numNewSub vector, newcb(1,i) is the
%          lower confidence bound for the i-th new subject, while
%          newcb(2,i) is the upper confidence bound.

function [newy newcb] = FPCpred(res, xx, yy, newx, new_tx, y, new_ty, FIT, alpha)

       if nargin < 9
           alpha = -1;
           if nargin < 8||isempty(FIT)||FIT == 0
               method = 'YPC';
           elseif FIT == -1
               method = 'YCO';
           end
       elseif isempty(alpha)
           alpha = .05;
       end               
       
       if isempty(yy)
           isYFun = 0;
       else isYFun = 1;
       end
       b = getVal(res, 'b');
       if isYFun == 1
            [K_x K_y] = size(b);
       else
           K_x = length(b);
           K_y = [];
       end
       [ypred, newxscore] =  FPCApred(xx, newx, new_tx);
       [newy newcb] = predict(newxscore, new_tx, xx, y, yy, new_ty, b, isYFun, K_x, K_y, method, alpha); 
end
