function [t,s,gradArea,plateauArea,plateauLimits] = grEncoding(fov,numSamples,...
    gradAmplitude,gradSlewrate,tstep,gamma)
%
% SEQUENCE.WAVEFORM.CARTESIANENCODING
%
%	Function that generates a cartesian encoding gradient,
%   based on FOV, resolution and system specs.
%   Area of the gradient is computed so that it covers 1/2 of K-space
%   (Intended to generate the Phase Encodings).
%
% INPUT
%   fov             field of view, in m
%   numSamples      number of samples to take during the plateau
%   gradAmplitude   amplitude of the gradient, T/m
%   gradSlewrate    slew rate of the gradient, T/m/s
%   tstep           time discretization
%   gamma           gyromagnetic ratio
%
% OUTPUT
%   t               discretized time vector, starts in 0
%   s               signal vector
%   gradArea        total area of the gradient
%   plateauArea     area in plateau region
%   plateauLimits   [start,end] points of the plateau
%
%========================  CORSMED AB © 2020 ==============================
%

if (nargin < 1 || isempty(fov))
    fov=0.256;
end
if (nargin < 2 || isempty(numSamples))
    numSamples=256;
end
if (nargin < 3 || isempty(gradAmplitude))
    gradAmplitude=0.030;
end
if (nargin < 4 || isempty(gradSlewrate))
    gradSlewrate=0.0;
end
if (nargin < 5 || isempty(tstep))
    tstep = 1e-7;
end
if (nargin < 6 || isempty(gamma))
    gamma = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
end

% encoding area (T*s) that we need for Kmax (half of full k-space)
%  (scaled in tstep temporal units)
gradArea = numSamples/(gamma*fov*2)/tstep;

% call function to generate trapezoidal based on area
[t,s,gradArea,plateauArea,plateauLimits] = ...
    simulatorBis.sequence.waveforms.grTrapArea(gradArea,gradAmplitude,gradSlewrate,tstep);
