function [obj_fun]=expS(par,info)
%y=par(1)*exp(-par(2)*info.t);
[fx, Gx] = fwdFun(par.x);

obj_fun=fx;
%obj_fun=(info.z-y)./info.sd;
