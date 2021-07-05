function [B1,Dph,Dfr] = addRF_AdiabaticTanTanh_FullPassage(B1max,...
    fmax,z,k,rf_duration,dt)
% It returns the rf amplitude B1(t) and phase Dph modulation functions
% of the adiabatic pulse as it is described in the paper by Kellman MRM, 
% 71:1428-1434 (2014) "Adiabatic Inversion Pulses for Myocardial T1 Mapping"
% It returns also the freq. modulation matrix which is zero-filled since
% this pulse is phase modulated only.

t       = 0:dt:rf_duration-dt;
t_hp    = 0:dt:rf_duration/2-dt;
t_rhp   = 0:dt:rf_duration/2-dt;

phmax = (-2*pi*fmax*rf_duration/(k*tan(k)))*log(cos(k));

% HALF PASSAGE 
% according to paper by Garwood, 1991 pg.515
rf_duration_hp  = rf_duration/2;  % ??
B1_hp           = B1max*tanh(z*t_hp/rf_duration_hp);
Dph_hp          = phmax-(2*pi*fmax*rf_duration_hp/(k*tan(k)))*...
    log(cos(k*(1-t_hp/rf_duration_hp))/cos(k));


% REVERSE HALF PASSAGE
% according to paper by Kellman, 2014, MRM 71:1428-1434
rf_duration_rhp = rf_duration/2;  % ??
B1_rhp          = B1max*tanh(z*(1-t_rhp)/rf_duration_rhp);
Dph_rhp         = phmax-(2*pi*fmax*rf_duration_rhp/(k*tan(k)))*log(cos(k*t_rhp/rf_duration_rhp)/cos(k));

B1  = [B1_hp,fliplr(B1_hp)];
Dph = [Dph_hp,Dph_rhp];
Dfr = zeros(size(Dph));

end