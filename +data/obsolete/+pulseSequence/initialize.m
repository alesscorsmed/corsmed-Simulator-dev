function [pulseSequence] = initialize()
%
% DATA.PULSESEQUENCE.INITIALIZE
%
%	Function that initializes a pulse sequence data structure.
%   Returns a pulseSequence with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   pulseSequence   pulseSequence structure
%
%========================  CORSMED AB © 2020 ==============================
%

%% basic info
pulseSequence.name       = 'EmptySequence';
pulseSequence.type       = 'N/A';
pulseSequence.endEvent   = 'relax'; % indicates what happens at the end
pulseSequence.totalTime  = 0.0; % total time in seconds
pulseSequence.numSteps   = 0;   % number of time steps
pulseSequence.numRxs     = 0;   % number of readout points
pulseSequence.gamma      = 42.577478518e6; 


%% main sequence data
% waveforms
pulseSequence.time      = zeros(pulseSequence.numSteps,1); % times
pulseSequence.timeDiff  = zeros(pulseSequence.numSteps,1); % time deltas
pulseSequence.rxSignal  = zeros(pulseSequence.numSteps,1); % receiver readout
pulseSequence.swcSignal = zeros(pulseSequence.numSteps,1); % software crusher
pulseSequence.gxSignal  = zeros(pulseSequence.numSteps,1); % x gradient
pulseSequence.gySignal  = zeros(pulseSequence.numSteps,1); % y gradient
pulseSequence.gzSignal  = zeros(pulseSequence.numSteps,1); % z gradient
pulseSequence.rfmSignal = zeros(pulseSequence.numSteps,1); % RF magnitude
pulseSequence.rfpSignal = zeros(pulseSequence.numSteps,1); % RF phase
pulseSequence.rffSignal = zeros(pulseSequence.numSteps,1); % RF frequency

% for area compressed sequences
pulseSequence.hArea     = zeros(pulseSequence.numSteps,1); % time deltas for areas
pulseSequence.gxArea    = zeros(pulseSequence.numSteps,1); % x gradient area delta
pulseSequence.gyArea    = zeros(pulseSequence.numSteps,1); % y gradient area delta
pulseSequence.gzArea    = zeros(pulseSequence.numSteps,1); % z gradient area delta

% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts     = 1; % number of parts
pulseSequence.partType{1}  = 'N/A'; % type of part: RF / RO / GR / DF
pulseSequence.partLimits   = zeros(pulseSequence.numParts,2); % index start/end of part

% for indicating RO parts
pulseSequence.rxLimits     = [];     
