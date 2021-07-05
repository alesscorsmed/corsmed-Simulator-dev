function [pulseSequence] = initializeSequence()
%
% DATA.SIMULATION.INITIALIZESEQUENCE
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
pulseSequence.type      = 'N/A';
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.gamma     = 42.577478518e6; 
pulseSequence.totalTime = 0.0; % total time in seconds
pulseSequence.numSteps  = 0; % number of time steps
pulseSequence.numRxs    = 0; % number of readout points
pulseSequence.numReps   = 0; % number of Total repetitions
pulseSequence.numEnc    = 0; % number of encoding repetitions
pulseSequence.numShots  = 0; % number of shots
pulseSequence.numTL     = 0; % number of reps in shot
pulseSequence.TE        = 0.0; % repetition TE;
pulseSequence.TR        = 0.0; % repetition TR;
pulseSequence.TI        = 0.0; % preparation TI;
pulseSequence.ESP       = 0.0; % echo spacing
pulseSequence.effTE     = 0.0; % effective echo time
pulseSequence.effTR     = 0.0; % effective TR
pulseSequence.effTI     = 0.0; % effective TI (to center of K space)
pulseSequence.maxIL     = 1; % maximum interleaving possible (slices in a TR)

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

% pulsed gradients for diffusion
pulseSequence.encPG     = []; % info for diffusion
pulseSequence.gdwSignal = []; % diffusion gradients

% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts     = 1; % number of parts
pulseSequence.partType{1}  = 'N/A'; % type of part: RF / RO / GR / DF
pulseSequence.partLimits   = zeros(pulseSequence.numParts,2); % index start/end of part

% for indicating RO parts
pulseSequence.rxLimits     = []; % stacked numRO x 2 w/ start and end indexes of RO samples   
