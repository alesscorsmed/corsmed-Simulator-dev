function [t,s,gradArea,plateauArea,plateauLimits] = grTrapPlateau(plateauTime,...
    gradAmplitude,gradSlewrate,tstep)
%
% SEQUENCE.WAVEFORM.GRTRAPEZOID
%
%	Function that generates a trapezoidal gradient,
%   based on plateau time and amplitude.
%
% INPUT
%   plateauTime     time in plateau, in s
%   gradAmplitude   amplitude of the gradient, T/m
%   gradSlewrate    slew rate of the gradient, T/m/s
%   tstep           time discretization, in s
%
% OUTPUT
%   t               discretized time vector, starts in tstep
%   s               signal vector
%   gradArea        total area of the gradient
%   plateauArea     area in plateau region
%   plateauLimits   [start,end] points of the plateau
%
%========================  CORSMED AB © 2020 ==============================
%

if (nargin < 1 || isempty(plateauTime))
    plateauTime=3e-3;
end
if (nargin < 2 || isempty(gradAmplitude))
    gradAmplitude=0.030;
end
if (nargin < 3 || isempty(gradSlewrate))
    gradSlewrate=150.0;
end
if (nargin < 4 || isempty(tstep))
    tstep = 1e-7;
end


% generate the time discretization and signal
duration = round(plateauTime/tstep);
if (gradSlewrate <= 0.0)
    % assume ideal waveforms, with no slewrate
    riseTime = 0;
else
    % take into account encoding area (T*s) covered by 
    % the rise and fall of the gradient (to max) due to SR   
    riseTime = ceil((gradAmplitude/gradSlewrate)/tstep);
    duration = duration + 2*riseTime;
end
t = 1:duration;
s = gradAmplitude*ones(size(t));
if riseTime > 0
    s(1:riseTime+1) = (t(1:riseTime+1)-1)*(gradAmplitude/riseTime);
    s(duration-riseTime:duration) = s(duration-riseTime:duration) - s(1:riseTime+1);
end
% numerical area
gradArea = sum(s)*tstep;
% numerical plateau area
plateauLimits = [riseTime+1, duration-riseTime];
plateauArea = sum(s(riseTime+1:duration-riseTime))*tstep;
% scale time by time steps
t = t*tstep;

