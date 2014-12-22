%function createDesignPlot(datafile, isColorPlot, noDiagonal, isIndicator, yname)
%This function creates the design plot of the data
%Input: 
%datafile:  i) character string of the file name for input time matrix,
%              where missing values are denoted by "NaN"
%           ii) input time cell array, same as "t" from PCA()
%isColorPlot: 1  create color plot
%             0  create black and white plot
%noDiagonal:  1  remove diagonal time pairs 
%             0  do not remove diagonal time pairs
%isIndicator: 1 each element of the count matrix is an indicator of a
%               presence of time pair
%             0 each element of the count matrix is the total number of 
%               presences of a given time pair
%yname:       same as "yname" in PCA()
function createDesignPlot(datafile, isColorPlot, noDiagonal, isIndicator, yname)
    
     if ischar(datafile) == 1
       timeMatrix = load(datafile); %time matrix, where missing is denoted by NaN
       t  = {};
       for i = 1:size(timeMatrix,1)
         tmp = timeMatrix(i,:);
         t{i} = tmp(~isnan(tmp));   %convert it to cell array
       end
     else
       t = datafile;                %no conversion is needed, if the input is cell array
     end
     out1 = unique(cell2mat(t));
     [res] = designPlotCount(t, out1, noDiagonal, isIndicator);
     if isColorPlot == 1
       createColorPlot(res,out1);
     else
       createBlackPlot(res,out1);
     end
     
     axis([min(out1) max(out1) min(out1) max(out1)]);
     xlabel('$T_{ij}$', 'Interpreter', 'latex');
     ylabel('$T_{ik}$', 'Interpreter', 'latex');
     title(['\bf{Design Plot' ' of ' yname '}']);
     hold off
end

function createBlackPlot(res, out1)
     figure;
     for i = 1:length(out1)

         idx = find(res(i,:) > 0);
         plot(out1(i)*ones(1, length(idx)), out1(idx),'k.');
         hold on;

     end
     
end

function createColorPlot(res, out1)
     figure;
     for i = 1:length(out1)
       tmp = res(i,:); 
       idx = find(tmp > 0);
       for j = 1:length(idx)
          plot(out1(i), out1(idx(j)), searchCol(tmp(idx(j))), 'Marker', '.', 'LineStyle', 'none');
          hold on;
       end

    end

end

function [col] = searchCol(val) 

   if val == 1
       col = 'r';
   elseif val == 2
       col = 'm';
   elseif val >= 3 && val <= 5
       col = 'g';
   elseif val > 6
       col = 'b';
   end


end
