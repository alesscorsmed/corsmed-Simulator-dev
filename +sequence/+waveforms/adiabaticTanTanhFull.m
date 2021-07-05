function [time,signalM,signalP,signalF,BW]= adiabaticTanTanhFull(...
    duration,b1Max,fMax,z,k,tstep,gamma)
%
% SEQUENCE.WAVEFORM.ADIABATICTANTANHFULL
%
%   It returns the rf amplitude B1(t) and phase Dph modulation functions
%   of the adiabatic pulse as it is described in the paper by Kellman MRM, 
%   71:1428-1434 (2014) "Adiabatic Inversion Pulses for Myocardial T1 Mapping"
%   It returns also the freq. modulation matrix which is zero-filled since
%   this pulse is phase modulated only.
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
if (nargin < 1|| isempty(duration))
    duration=3e-3;
end
if (nargin < 2 || isempty(b1Max))
    b1Max = 0.000015;
end
if (nargin < 3 || isempty(fMax))
    fMax = 9500;
end
if (nargin < 4 || isempty(z))
    z = 10;
end
if (nargin < 5 || isempty(k))
    k = 1.525373047373320;
end
if (nargin < 6 || isempty(tstep))
    tstep = 1e-6;
end
if (nargin < 7 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end

%% create basic set up and time discretization
totalTimesteps      = 2*round(duration/2/tstep);
time                = tstep*(1:totalTimesteps);
signalM             = 0*time;
signalP             = 0*time;
signalF             = 0*time;
BW                  = 0;

% parameters
halfTimesteps       = round(duration/2/tstep);
durationHalfPass    = tstep*halfTimesteps;
timeHalfPass        = tstep*(1:halfTimesteps);
phMax               = -2*pi*fMax*duration/(k*tan(k)) *log(cos(k));

%% half passage 
% according to paper by Garwood, 1991 pg.515
b1HalfPass          = b1Max*tanh(z*timeHalfPass/durationHalfPass);
phHalfPass          = phMax ...
    - (2*pi*fMax*durationHalfPass/(k*tan(k)))*...
    log(cos(k*(1-timeHalfPass/durationHalfPass))/cos(k));

%% reverse half passage
% according to paper by Kellman, 2014, MRM 71:1428-1434
b1ReverseHalfPass	= b1Max*tanh(z*(1-timeHalfPass/durationHalfPass));
phReverseHalfPass 	= phMax ...
    -(2*pi*fMax*durationHalfPass/(k*tan(k)))* ...
    log(cos(k*timeHalfPass/durationHalfPass)/cos(k));

%% combine both
signalM(:)  = [ b1HalfPass(:); b1ReverseHalfPass(:) ];
signalP(:)  = [ phHalfPass(:); phReverseHalfPass(:) ] ;

