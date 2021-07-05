function [t,s,gradArea,plateauArea,plateauLimits] = grTrapArea(gradArea,...
    gradAmplitude,gradSlewrate,tstep)
%
% SEQUENCE.WAVEFORM.GRTRAPAREA
%
%	Function that generates a trapezoidal gradient,
%   based on gradient area.
%
% INPUT
%   gradArea        total area of the gradient
%   gradAmplitude   amplitude of the gradient, T/m
%   gradSlewrate    slew rate of the gradient, T/m/s
%   tstep           time discretization, in s
%
% OUTPUT
%   t               discretized time vector, starts in tstep
%   s               signal vector
%   gradArea        total area of the gradient (numerically corrected)
%   plateauArea     area in plateau region
%   plateauLimits   [start,end] points of the plateau
%
%========================  CORSMED AB © 2020 ==============================
%

if (nargin < 1 || isempty(gradArea))
    gradArea=3e-3*0.030;
end
if (nargin < 2 || isempty(gradAmplitude))
    gradAmplitude=0.030;
end
if (nargin < 3 || isempty(gradSlewrate))
    gradSlewrate=0.0;
end
if (nargin < 4 || isempty(tstep))
    tstep = 1e-7;
end

if (gradSlewrate <= 0.0)
    % assume ideal waveforms, with no slewrate
    riseTime = 0;
else
    % take into account encoding area (T*s) covered by 
    % the rise and fall of the gradient (to max) due to SR   
    riseTime = ceil((gradAmplitude/gradSlewrate)/tstep);
end
riseArea = gradAmplitude*riseTime*tstep;
% area and time we need to spend in plateau
plateauArea = max(gradArea-riseArea,0);
plateauTime = ceil(plateauArea/gradAmplitude/tstep);

% generate the time discretization and signal
duration = plateauTime + 2*riseTime;
t = 1:duration;
s = gradAmplitude*ones(size(t));
% if we have slew rate, create the ramps
if riseTime > 0
    s(1:riseTime+1) = (t(1:riseTime+1)-1)*(gradAmplitude/riseTime);
    s(duration-riseTime:duration) = s(duration-riseTime:duration) - s(1:riseTime+1);
end

% scale to correct any numerical area approximations
s = s*gradArea/(sum(s)*tstep);
gradArea = sum(s)*tstep;
% numerical plateau area
plateauLimits = [riseTime+1, duration-riseTime];
plateauArea = sum(s(riseTime+1:duration-riseTime));
plateauArea = plateauArea*tstep;
% scale time by time steps
t = t*tstep;
