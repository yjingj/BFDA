% Predict the response value given the FPC scores of a new predictor
% function.
function fitted_y=predict(xi,coeff)
n=size(xi,1);
qdm1=[ones(n,1),qdmscores1(xi)];
alpha=getVal(coeff,'alpha');
beta=getVal(coeff,'beta');
gamma=getVal(coeff,'gamma');
fitted_y=[alpha',beta',gamma']*qdm1';
