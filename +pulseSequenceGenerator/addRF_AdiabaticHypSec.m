function [B1,Dph,Dfr] = addRF_AdiabaticHypSec(B1max,Dwmax,rf_duration,dt)
% It returns the rf amplitude B1(t) and frequency Dfr(t) modulation functions
% of the adiabatic pulse as it is described in the paper with the title
% "Adiabatic RF pulses: in vivo applications", pg 250-252

t       = 0:dt:rf_duration-dt;

beta    = asech(0.01);

B1      = B1max*sech(beta*((2*t/rf_duration)-1));
Dfr     = Dwmax*tanh(beta*(1-(2*t/rf_duration)));
Dph     = zeros(size(Dfr));

end