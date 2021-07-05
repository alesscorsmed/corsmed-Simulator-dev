function [F,J] = fitFun2(x, xdata)
% Fitting function w.r.t apparent T1 == x(2)
myExp = exp(-xdata./x(2));
Fp = x(1) - x(3) .* myExp;
F = abs(Fp); 

if nargout > 1   
   J = [ones(size(xdata));  -x(3) .* myExp .* xdata ./ x(2)^2 ; - myExp ]' .* sign(Fp)';
end