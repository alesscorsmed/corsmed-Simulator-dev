function [time,signalM,signalP,signalF,BW]=rfGauss(...
    duration,angle,widthFactor,tstep,gamma)
%
% SEQUENCE.WAVEFORM.RFGAUSS
%
%	Function that generates a Gaussian RF pulse.
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
if (nargin < 3 || isempty(widthFactor))
    widthFactor=0.15; % approximated "width" of the RF pulse w.r.t duration
end
if (nargin < 4 || isempty(tstep))
    tstep = 1e-6;
end
if (nargin < 5 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end

%% create basic set up and time discretization
totalTimesteps	= round(duration/tstep);
duration     	= tstep*totalTimesteps;
time        	= tstep*(1:totalTimesteps);
signalP     	= 0*time;
signalF     	= 0*time;
%% compute normalized Gaussian signal
mu      = duration/2;
tau 	= widthFactor*duration;
signalM = 1/(tau*sqrt(2*pi)) * exp(-(time - mu).^2 /(2*tau^2));
%% scale RF magnitude for flip angle
area    = sum(signalM);
B1magn  = angle/(360*gamma*area*tstep);
signalM = B1magn*signalM;
%% BW
fwhm = 2*sqrt(2*log(2))*tau; % full-width half magnitude
BW   = 0.44/fwhm; % approx
