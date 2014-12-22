%function [newy, newt] = binData(y,t,regular, verbose, numBins)
%This function is to prebin the data.
%
% ======
% Input:
% ====== 
%      y : 1 * n cell array of repeated measurements
%          where y{i} is a row vector of observed 
%          repeated measurements for subject i.
%
%          Same one as in FPCA.m or PCA.m
%
%      t : 1 * n cell array of time points
%          where t{i} is a row vector of observed
%          time points for subject i
%
% regular : 0  sparse and irregular data
%           1  regular data with missing values
%           2  completely balanced data
%
% verbose: 'on' display diagnostic messages       
%          'off' suppress diagnostic messages
%
% numBins: i) positive integer denotes the number of bins define 
%          in the time domain of all t{i}'s
%          ii) if numBins is not specified:
%                 
%            --  regular = 1 or 2
%                m = max of n_i, where n_i is the number of repeated measurements
%                for ith subject.
%            --  regular = 0
%                m = median of n_i
%
%                When n <= 5000 subjects and m <= 400 or m <= 20, no binning.
%                When n <= 5000 subjects and m > 400, numBins = 400.
%                When n > 5000 subjects, compute
%
%                m* = max(20, (5000-n)*19/2250+400);
%
%                if m > m*, numBins = m*
%                if m <= m*, no binning
%
%
% ======
% Output:
% ======
%   newy : 1 * n cell array of repeated measurements after binning
%   newt : 1 * n cell array of corresponding time points, t{i} are 
%          the midpoints of each bin
%
%   newy and newt will be [] if the median number of repeated
%   measurements, i.e., median of n_i <= 20.
%  
%    ex: assume y and t have been defined
%      1) Use default number of bins
%        regular = 0;   %here, regular can be either 0,1 or 2
%        [newy, newt] = binData(y,t,regular,'off')
%    or 
%      2) Use 10 bins
%        [newy, newt] = binData(y,t,[],'off',10)       
%
function [newy, newt] = binData(y,t,regular, verbose, numBins)
    
    n = length(t);
    ni = zeros(1,n);
    for i = 1:n
         ni(i) = length(t{i});
    end

    if isempty(regular) 
       regular = 0;
    end

    if isempty(verbose)
       verbose = 'on';
    end

    if regular == 0
      m = median(ni);
    else
      m = max(ni);
    end

    if nargin == 4 
        numBins = getBINnum(n,m,regular, verbose);
    elseif numBins <= 0
        fprintf(1,'Warning: number of bins must be positive integer! Reset to default number of bins now!\n');
        numBins = getBINnum(n,m,regular, verbose);
    end
    
    if isempty(numBins)
       newy = []; newt = [];
       return;
    end
    
    numBins = ceil(numBins);
    
    newt = cell(1,n);
    newy = cell(1,n);
    tt = cell2mat(t);  
    a0 = min(min(tt));
    b0 = max(max(tt));

    if strcmp(verbose, 'on') == 1
       fprintf(1,'Start binning: \n');
    end

    for i = 1:n
      res = binning(t{i}, y{i},numBins, 1,1,a0,b0);
      newt{i} = res.midpoint;   %or getVal(res, 'midpoint');
      newy{i} = res.newy;       %or getVal(res, 'newy');
    end
    if strcmp(verbose, 'on') == 1
       fprintf(1,['Number of bins selected for the input data in the time domain [' num2str(a0) ',' num2str(b0) '] is ' num2str(numBins) ' .\n']);
    end
end

function numBin = getBINnum(n,m,regular,verbose)

   numBin = [];
   if m <= 20
     if regular == 0
        str = 'Median of ni';
     else
        str = 'Maximum of ni';
     end
     if strcmp(verbose, 'on') == 1
        fprintf(1, [str ' is no more than 20!No binning is performed!\n']);
     end
     return;
   end
  
   if m > 400
     numBin = 400;
   end

   if n > 5000
      mstar = max(20,(5000-n)*19/2250+400);
      if m > mstar
          numBin = ceil(mstar);
      else
          if strcmp(verbose, 'on') == 1
             fprintf(1,'No binning is needed.\n');
          end
          numBin = [];
          return;
      end
   end
   
   if isempty(numBin) && strcmp(verbose, 'on') == 1
      fprintf(1,'No binning is needed.\n');
   end
 
end
