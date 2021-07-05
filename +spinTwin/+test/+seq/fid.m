function [pulseSequence] = fid(totalTime, dt )
%
% SPINTWIN.TEST.SEQ.FID
%
%	Generates a "FID sequence"
%
% INPUT
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:test:seq:fid';
% args
if (nargin<2) || isempty(dt)
    totalTime = 1e-3; 
end
if (nargin<2) || isempty(dt)
    dt = 1e-4; 
end

% constant
gamma = 42.577478518e6;

numSteps = ceil(totalTime/dt);

%% empty sequence
pulseSequence = data.simulation.initializeSequence();
pulseSequence.type      = 'FID';
% waveforms
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.timeDiff  = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1); % receiver readout
pulseSequence.swcSignal = zeros(numSteps,1); % software crusher
pulseSequence.gxSignal  = zeros(numSteps,1); % x gradient
pulseSequence.gySignal  = zeros(numSteps,1); % y gradient
pulseSequence.gzSignal  = zeros(numSteps,1); % z gradient
pulseSequence.rfmSignal = zeros(numSteps,1); % RF magnitude
pulseSequence.rfpSignal = zeros(numSteps,1); % RF phase
pulseSequence.rffSignal = zeros(numSteps,1); % RF frequency
pulseSequence.gdwSignal = zeros(numSteps,3); % diffusion gradients
% assign blocks
pulseSequence.time(:)       = dt*(1:numSteps);
pulseSequence.timeDiff(:)   = dt;
pulseSequence.gzSignal(:)   = 0.0;
pulseSequence.rxSignal(:)   = 1:numSteps;
pulseSequence.totalTime     = dt*numSteps; % total time in seconds
pulseSequence.numSteps      = numSteps;   % number of time steps
pulseSequence.numRxs        = numSteps;   % number of readout points
% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts      = 1; % number of parts
pulseSequence.partType{1}   = 'RO'; % type of part: RF / RO / GR / DF
pulseSequence.partLimits    = [1, numSteps]; % index start/end of part
% for indicating RO parts
pulseSequence.rxLimits      = [1, numSteps];
