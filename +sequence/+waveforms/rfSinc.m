function [time,signalM,signalP,signalF,BW]=rfSinc(...
    duration,angle,cycles,tstep,gamma)
%
% SEQUENCE.WAVEFORM.RFSINC
%
%	Function that generates a sinc RF pulse.
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
if (nargin < 2 || isempty(angle))
    angle=90;
end
if (nargin < 3 || isempty(cycles))
    cycles=2;
end
if (nargin < 4 || isempty(tstep))
    tstep = 1e-6;
end
if (nargin < 5 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end

%% create basic set up and discretization
sideLobes      = (cycles-1)*2;
totalLobes     = sideLobes + 2;  % because the main lobe has width of 2
                                 % times the width of the sidelobe
totalTimesteps = round(duration/tstep);
sincStep       = totalLobes/totalTimesteps; % because the sinc function is
                                            % created in cycles. So, since I have 4
                                            % cycles in total, I divide this by the
                                            % timesteps that I want my sinc function to
                                            % have

% time and signal
tsinc   = 0:sincStep:totalLobes-sincStep;
signalM = sinc(tsinc-(2+sideLobes)/2);
time    = tstep*(1:totalTimesteps);

%% scale RF magnitude for flip angle
area    = sum(signalM);
B1magn  = angle/(360*gamma*area*tstep);
signalM = B1magn.*signalM; % signal magnitude
signalP = 0*signalM;       % signal phase shift
signalF = 0*signalM;       % signal frequency shift

%% Accumulated angle of the RF pulse
anglePerTimestep  = 360*signalM*gamma*tstep;
accumAngle        = cumsum(anglePerTimestep);

%% Calculate the RF BW
% More info on how to calculate the bandwidth (BW) of the rf pulse (sinc)
% see pg.37 on the book "Handbook of MRI Pulse Sequences"
t0 = duration/(cycles*2);
BW = 1/t0;
