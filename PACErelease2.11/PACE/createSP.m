%function createSP(FVE, no_opt, yname)
%This function creates the screePlot based on the
%FVE and no_opt from FPCA().
%Input FVE:  1*ngrid vector of fraction of variance explained
%Input no_opt: integer, automatically or subjectively selected K (number of
%              FPC scores and eigenfunctions.
%Input yname:  a character string of the name of the current function 
%              to be estimated. If missing, then it is set to be 'y'
function createSP(FVE, no_opt, yname)
if nargin < 3
  yname = 'y';
end

figure;
plot(0:length(FVE), [0 FVE]*100, '--ro', 'LineWidth',2,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
xlabel('\bf{No. of Principal Components}');
ylabel('\bf{FVE (%)}');
title(['\bf{Fraction of variance explained by No. of PC for function ' yname '}'])
hold on;
plot(linspace(0,no_opt,30), ones(30,1)*FVE(no_opt)*100,'b',...
     ones(30,1)*no_opt, linspace(0,FVE(no_opt)*100,30),'b');
text(no_opt+0.2, FVE(no_opt)*100-10, {[' k = ' num2str(no_opt) ', FVE = ' num2str(roundoff(FVE(no_opt)*100,3)) '%'] ' (final choice)'});
axis([0 length(FVE)+1 0 101]);
hold off;

end
