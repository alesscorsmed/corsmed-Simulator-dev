function [time,signalM,signalP,signalF,BW]= adiabaticBIR4(...
    duration,angle,b1Max,dwMax,beta,kappa,tstep,gamma)
%
% SEQUENCE.WAVEFORM.ADIABATICBIR4
%
%   According to Inv.Rad 1990 25:559-567 Staewen et al.
%
%   B1max and Dwmax for a 5ms IR pulse based on the files AM_BIR4.txt and
%   FM_BIR4.txt at your personal documents. It doesn't work for other IR
%   durations.
%
% INPUT
%   duration        total RF duration, in s
%   cycles          number of cycles of the sinc
%   angle           Flip angle, in degrees
%   tstep           time discretization
%   gamma           gyromagnetic ratio
%
% OUTPUT
%   time            discretized time vector, starts in tstep
%   signalM         magnitude signal vector
%   signalP         phase shift signal vector
%   signalF         frequency shift signal vector
%   BW              RF bandwidth, in Hz
%
%========================  CORSMED AB © 2020 ==============================
%
if (nargin < 1 || isempty(duration))
    duration=5e-3;
end
if (nargin < 2 || isempty(angle))
    angle = 180;
end
if (nargin < 3 || isempty(b1Max))
    b1Max = 0.000015;
end
if (nargin < 4 || isempty(dwMax))
    dwMax = 8489;
end
if (nargin < 5 || isempty(beta))
    beta = 5;
end
if (nargin < 6 || isempty(kappa))
    kappa = atan(10); %1.4711;
end
if (nargin < 7 || isempty(tstep))
    tstep = 1e-6;
end
if (nargin < 8 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end

%% create basic set up and time discretization
totalTimesteps  = round(duration/tstep);
time            = tstep*(1:totalTimesteps);
duration        = tstep*totalTimesteps;
signalM         = 0*time;
signalF         = 0*time;
signalP         = 0*time;
BW              = 0;

%% compute the signals
phaseShift = 180+angle/2; 
% tau as an array
tau = time/(duration/4);
% tau < 1
idx = tau < 1;
signalM(idx) = tanh(beta*(1-tau(idx)));
signalF(idx) = tan(kappa*tau(idx))/tan(kappa);
% tau = 1
idx = tau == 1;
signalM(idx) = 0;
signalF(idx) = phaseShift/(360*tstep);
% tau>1 && tau<=2
idx = (tau > 1) & (tau <= 2);
signalM(idx) = tanh(beta*(tau(idx)-1));
signalF(idx) = tan(kappa*(tau(idx)-2))/tan(kappa); %+(phaseShift*pi/(180*dt));
% tau>2 && tau<3
idx = (tau > 2) & (tau < 3);
signalM(idx) = tanh(beta*(3 - tau(idx)));
signalF(idx) = tan(kappa*(tau(idx)-2))/tan(kappa); %+(phaseShift*pi/(180*dt));
% tau = 3
idx = tau == 3;
signalM(idx) = 0;
signalF(idx) = -phaseShift/(360*tstep);
% tau>3 && tau<4
idx = (tau > 3) & (tau < 4);
signalM(idx) = tanh(beta*(tau(idx)-3));
signalF(idx) = tan(kappa*(tau(idx)-4))/tan(kappa);

%% scaling
signalM = b1Max*signalM;
signalF = dwMax*signalF;

