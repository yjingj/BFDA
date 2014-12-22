%function [res] = getRawCov(y,t,out1new, mu,regular, error)
%obtain raw covariance
%Input y :       1*n cell array of the observed repeated measurements from n subjects
%Input t :       1*n cell array of the observed time points from n subjects
%Input out1new:  1*m vector of time points correspond to mu
%Input mu:       1*m vector of fitted mean functions from Step I, corresponding to
%                pooled unique time points from t
%Input regular:  same as in PCA.m
%Input error:    1 with measurement error assumption
%                0 without measurement error assumption
%
%Output res: a struct object that contains tpairn, cxxn, indx,win and cyy
%    tpairn:  2 * N  vector denotes the  pairs of time points for subject 
%                concatenating as two vectors 
%                if error = 1, all (t_ij, t_ij) will be removed  
%      cxxn:    1 * N vector of raw covariance corresponding to tpairn      
%      indx:    1 * N vector of indices for each subject
%       win:    1 * N weight matrix for the 2-D smoother for covariance function
%       cyy:    1 * M vector of raw covariance corresponding to all pairs of time points,
%               i.e., it is the same as cxxn if error = 0
function [res] = getRawCov(y,t,out1new, mu,regular, error)
xx1 = []; xx2 = [];yy1 = []; yy2 = [];
indx = []; count = [];
ncohort = length(y);
out1  = unique(cell2mat(t));
mu = mapX1d(out1new, mu, out1);

  if regular == 2                                       %regular balance data

      yy = reshape(cell2mat(y),length(y{1}), ncohort)';
      MU = repmat(mu,ncohort,1); 
      t1 = t{1};
      clear t y;
      yy = yy-MU;
      cyy = yy'*yy/ncohort;
      cyy = cyy(:)';
      cxxn = cyy;
      [xx,yy] = meshgrid(t1);
      clear t1;
      tpairn = [xx(:)'; yy(:)'];
      clear xx yy;
      if error==1
        tneq=find(tpairn(1,:)~=tpairn(2,:));
        cxxn=cyy(tneq);
        tpairn=tpairn(:,tneq);
      end
      win = ones(1,length(cxxn));
 
 elseif regular == 1

      yy = zeros(ncohort, length(out1));
      ID = zeros(ncohort,length(out1));
      for i = 1:ncohort
         id = ismember(out1,t{i});
         yy(i,id == 1) = y{i};
         ID(i,:) = id;
      end
    
      MU = repmat(mu,ncohort,1);
      MU(ID == 0) = 0;
      yy = yy-MU;
      count = ID'*ID;
      %cyy = (yy'*yy)./count;
      
      cyy = yy'*yy;
      cyy = cyy(:)';
      count = count(:); 
      cyy = cyy(count ~= 0);
      cxxn = cyy; 

      [xx,yy] = meshgrid(out1);
      clear MU out1 t;
      tpairn = [xx(:)'; yy(:)'];
      tpairn = tpairn(:,count ~=0);
      clear xx yy;
      count = count(count ~= 0);
      if error==1
        tneq=find(tpairn(1,:)~=tpairn(2,:));
        cxxn=cyy(tneq);
        tpairn=tpairn(:,tneq);
        count = count(tneq);
      end
      win = ones(1,length(cxxn));


 else                                          %sparse irregular or regular with missing data
      for i = 1:ncohort

        [x1,x2] = meshgrid(t{i});
        [y1,y2] = meshgrid(y{i});
        xx1 = [xx1; x1(:)];                    %vectorize observed time points
        xx2 = [xx2; x2(:)];                    %vectorize observed time points
        indx = [indx i*ones(1,length(x1(:)))]; %index for subject i
        yy1 = [yy1; y1(:)];                    %vectorize observed measurements
        yy2 = [yy2; y2(:)];                    %vectorize observed measurements

      end
      [outgrid1, ignore, id1] = unique(xx1);
      [outgrid2, ignore, id2] = unique(xx2);
 
      cyy = (yy1'-mu(id1)).*(yy2'-mu(id2));
      tpairn = [xx1 xx2]';
      
      if error == 1
        tneq = find(xx1 ~= xx2);
        indx = indx(tneq);
        tpairn = tpairn(:,tneq);
        cxxn = cyy(tneq);     
      elseif error == 0
        cxxn = cyy;     
      end
      win = ones(1, length(cxxn));
    
 end
     res = struct('tpairn',tpairn, 'cxxn',cxxn, 'indx',indx, 'win',win,'cyy',cyy,'count',count);

end
