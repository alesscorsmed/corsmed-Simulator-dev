function [pulseSequence, grData] = pgReversed( grData, mrSystem, dt )
%
% SPINTWIN.TEST.SEQ.PGREVERSED
%
%	Generates a pair of reversed gradients (+G -G) separated by an interval
%
% INPUT
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:test:seq:pgReversed';
% args
if (nargin<1) || isempty(grData)
    grData.echoTime         = 10e-3;
    grData.duration      	= 5e-3;
    grData.amplitude     	= 0.040;
    grData.direction        = 'x'; % encoding direction 'x' / 'y' / 'z'
end
if (nargin<2) || isempty(mrSystem)
    mrSystem.b0             = 1.5000;
    mrSystem.maxGStrenght   = 0.0400;
    mrSystem.SlewRate       = 0.0;
end
if (nargin<3) || isempty(dt)
    dt = 1e-6;
end

% constant
gamma = 42.577478518e6;

%% generate Pulsed Gradient
[pgTime,pgSignal,pgArea,~,pgPlatLimits] = sequence.waveforms.grTrapPlateau(...
        grData.duration, grData.amplitude, mrSystem.SlewRate, dt);
% compute remaining data    
grData.area     = pgArea; % gradient area
grData.time     = pgTime(end); % total time
grData.plateau  = diff(pgTime(pgPlatLimits)); % plateau time
grData.beta     = (2*pi*gamma)^2 *(grData.area^2) *(grData.echoTime-grData.time/3);

    
%% combine the gradients into a sequence
tWait        = 10e-6;
numTimesGrad = length(pgTime);
numTimesWait = round(tWait/dt);
numTimesTAU  = max(round(grData.echoTime/dt), numTimesGrad+numTimesWait);

% start/end of forward gradient
iniFwGr         = numTimesWait;
endFwGr         = iniFwGr+numTimesGrad-1;
% start/end of rewind gradient
iniRwGr         = iniFwGr+numTimesTAU;
endRwGr         = iniRwGr+numTimesGrad-1;
% total time of sequence
numSteps        = ceil(iniRwGr + numTimesGrad/2 + numTimesTAU);

% update and compute remaining PG data
grData.echoTime = dt*numTimesTAU;


%% pulse sequence
pulseSequence = data.simulation.initializeSequence();
pulseSequence.type         = 'PG Reversed';
pulseSequence.gamma        = gamma;

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

%% assign blocks
pulseSequence.time(:)     = dt*(1:numSteps);
pulseSequence.timeDiff(:) = dt;
%% PG gradients (Fw and Rw)
switch lower(grData.direction)
    case 'x'
        pulseSequence.gxSignal(iniFwGr:endFwGr,1) =  pgSignal(:);
        pulseSequence.gxSignal(iniRwGr:endRwGr,1) = -pgSignal(:);
    case 'y'
        pulseSequence.gySignal(iniFwGr:endFwGr,1) =  pgSignal(:);
        pulseSequence.gySignal(iniRwGr:endRwGr,1) = -pgSignal(:);
    case 'z'
        pulseSequence.gzSignal(iniFwGr:endFwGr,1) =  pgSignal(:);
        pulseSequence.gzSignal(iniRwGr:endRwGr,1) = -pgSignal(:);
    otherwise
        msg = sprintf('gradient direction %s not supported',...
            grData.direction );
        ME = MException(['error:',functionName], '%s', msg);
        throw(ME);
end
pulseSequence.rxSignal(:) = 1:numSteps;
rxLimits = [1, numSteps];

%% pulse sequences info
pulseSequence.totalTime    = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps     = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs       = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.effTE        = grData.echoTime;
pulseSequence.rxLimits     = rxLimits; % starting / end indexes of Readouts

%% for splitting into parts (RF vs Non-RF)
numParts = 1;
partLimits = [1, numSteps];
partType{1} = 'RO';
pulseSequence.numParts   = numParts; % number of parts
pulseSequence.partType   = partType; % type of part
pulseSequence.partLimits = partLimits; % index start/end of part




