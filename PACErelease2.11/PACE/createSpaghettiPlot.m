%  function createSpaghettiPlot(y,t)
%    This function creates a ``Spaghetti'' plot from functional data that are given in cell
%    array form. 
%  
%   Input 
%   y:          1*n cell array, y{i} is the vector of measurements for the ith
%              subject, i=1,...,n.
%   t:          1*n cell array, t{i} is the vector of time points for the ith
%              subject, at which corresponding measurements y{i} are taken,
%              i=1,...,n.
%   
%   Example:
%   createSpaghettiPlot(y,t);   %Plot the original data, where y contains the responses and t contains the corresponding time points
%   yy = FPCA(y,t);
%   y_predOrig = getVal(yy, 'y_predOrig');
%   createSpaghettiPlot(y_predOrig,t);    %Plot the predicted functions as obtained from the PACE fitting procedure
%  
%   % Plot the first 20 subjects from the original data and the fitted functions for direct comparison
%   subplot(1,2,1);
%   createSpaghettiPlot(y(1:20),t(1:20));
%   subplot(1,2,2); 
%   createSpaghettiPlot(y_predOrig(1:20),t(1:20));
%  
%     See also PCA

   
function createSpaghettiPlot(y,t)

  cm = colormap;
  n = length(y);
  plot(t{1},y{1},'');
  hold on;
  for i =1:n
      [sortti,idx] = sort(t{i});
      yi = y{i};
      line(sortti,yi(idx),'Color',cm(mod(i,8)*8+4,:));
  end
  title('Spagetti Plot');
end