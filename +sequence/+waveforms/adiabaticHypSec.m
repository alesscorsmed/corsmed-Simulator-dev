function [time,signalM,signalP,signalF,BW]= adiabaticHypSec(...
    duration,b1Max,dwMax,tstep,gamma)
%
% XRONOS.WAVEFORM.ADIABATICHYPSEC
%
%   It returns the rf amplitude B1(t) and frequency Dfr(t) modulation functions
%   of the adiabatic pulse as it is described in the paper with the title
%   "Adiabatic RF pulses: in vivo applications", pg 250-252
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
    duration=3e-3;
end
if (nargin < 2 || isempty(b1Max))
    b1Max = 0.0000270606;
end
if (nargin < 3 || isempty(dwMax))
    dwMax = 1343.39;
end
if (nargin < 4 || isempty(tstep))
    tstep = 1e-6;
end
if (nargin < 5 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end

%% create basic set up and time discretization
totalTimesteps  = round(duration/tstep);
time            = tstep*(1:totalTimesteps);

%% parameters
beta    = asech(0.01);
tarray  = 2*time/duration;

%% generate the signals: magnitude, frequency shift, phase shift
signalM     = b1Max*sech( beta*(tarray-1) );
signalF     = dwMax*tanh( beta*(1-tarray) );
signalP     = 0*time;

%% Calculate the RF BW
BW = 2*dwMax;
