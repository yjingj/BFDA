function [res invalid]=getSC(x, t_x, p_x, y, t_y, p_y, bwccov, nsvd, regular)

if isempty(p_x)
    p_x = setOptions('regular', regular, 'selection_k', nsvd, 'screePlot',0);
else
    p_x = setVal(p_x, 'selection_k', nsvd, 'screePlot',0);
end
xx = FPCA(x,t_x,p_x);
if isempty(p_y)
    p_y = setOptions('regular', regular, 'selection_k', nsvd, 'screePlot',0);
else
    p_y = setVal(p_y, 'selection_k', nsvd, 'screePlot',0);
end
yy = FPCA(y,t_y,p_y);

[res invalid] = getSC1(x, t_x, xx, y, t_y, yy, bwccov, nsvd, regular);

end