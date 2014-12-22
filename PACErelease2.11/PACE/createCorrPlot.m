%function createCorrPlot(yy)
%This function creates the correlation surface plot based on the
%results from FPCA() or FPCder()
%Input yy : returned object from FPCA().
%example:
%yy = FPCA(y,t,p);
%createCorrPlot(yy)
%or
%p = setDerOptions('nder',0:2);
%yy = FPCder(y,t,p);
%createCorrPlot(yy)  
function createCorrPlot(yy)

  no_opt = getVal(yy,'no_opt');
  if isempty(no_opt)
      fprintf(1,'Warning: Input is not a valid FPCA or FPCder output.');
  elseif no_opt == 1
      fprintf(1,'Warning: Correlation surface is not available when only one principal component is used.');
  else
      ops = getVal(yy,'ops');
      yname = getVal(ops,'yname');
      out21 = getVal(yy,'out21');
      xcorr = getVal(yy,'xcorr');
      
      figure;
      mesh(out21,out21,xcorr);
      xlabel('\bf{t}');
      ylabel('\bf{t}');
      zlabel('\bf{Correlation}');
      title(['\bf{Fitted correlation surface for function ' yname '}']);
  end
end
