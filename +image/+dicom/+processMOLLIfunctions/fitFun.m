function [F,J] = fitFun(x, xdata)
% Fitting function % w.r.t T1 == x(2) .* (x(3)./x(1)-1))
myExp = exp(-xdata./x(2) .* (x(3)./x(1)-1));
Fp = x(1) - x(3) .* myExp;
F = abs(Fp);  

if nargout > 1   
    Jx1 = ones(size(xdata)) - myExp .* xdata .* x(3).^2 ./ x(1)^2 ./ x(2);
    Jx2 = -x(3) .* myExp .* xdata ./ x(2)^2 .* (x(3)./x(1)-1) ;
    Jx3 = myExp .* (-1 + xdata .* x(3) ./ x(1) ./ x(2));
    J = [Jx1;  Jx2; Jx3]' .* sign(Fp)';
end