function [pulseSequence] = singleRF( rfData, mrSystem, dt )
%
% SPINTWIN.TEST.SEQ.SINGLERF
%
%	Generates a RFBlock with optional slice selection,
%   and pre/post gradient rewinds.
%
% INPUT
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:test:seq:singleRF';
% args
if (nargin<1) || isempty(rfData)
    rfData.type             = 'sinc';
    rfData.flipAngle      	= 90; % in degrees
    rfData.phase         	= 0.0; % in rads
    rfData.duration         = 1e-3;
    rfData.cycles           = 2;
    rfData.sliceThickness   = 3e-3;
    rfData.doSliceSelect 	= 0;
    rfData.preRWScale    	= 0;
    rfData.postRWScale    	= 0;
    rfData.keepTimingSS     = 0;
end
if (nargin<2) || isempty(mrSystem)
    mrSystem.b0             = 1.5000;
    mrSystem.maxGStrenght   = 0.0400;
    mrSystem.SlewRate       = 0.0;
end
if (nargin<3) || isempty(dt)
    dt=1e-6;
end

%% generate the waveforms for the RF
gamma = 42.577478518e6;
[time,rfm,rfp,rff,gr,rfLimits] = sequence.excitations.sliceSelectRF(...
    rfData, mrSystem.maxGStrenght, mrSystem.SlewRate, gamma, dt );

%% sequences data
pulseSequence = data.simulation.initializeSequence();
pulseSequence.type         = 'Single RF';
pulseSequence.gamma        = gamma;
pulseSequence.totalTime    = time(end); % total time in seconds
pulseSequence.numSteps     = length(time); % number of time steps
pulseSequence.numRxs       = 1; % number of readout points

% waveforms
pulseSequence.time      = time;
pulseSequence.timeDiff  = diff([0; time]);
pulseSequence.rxSignal  = zeros(pulseSequence.numSteps,1); % receiver readout
pulseSequence.rxSignal(end) = 1; % unique readout at end
pulseSequence.swcSignal = zeros(pulseSequence.numSteps,1); % software crusher
pulseSequence.gxSignal  = zeros(pulseSequence.numSteps,1); % x gradient
pulseSequence.gySignal  = zeros(pulseSequence.numSteps,1); % y gradient
pulseSequence.gzSignal  = gr; % z gradient
pulseSequence.rfmSignal = rfm; % RF magnitude
pulseSequence.rfpSignal = rfp; % RF phase
pulseSequence.rffSignal = rff; % RF frequency

% for splitting into parts (RF vs Non-RF)
numParts = 1;
if rfLimits(1) > 1
    numParts = numParts + 1;
    partLimits = [1, rfLimits(1)-1; rfLimits];
    partType{1} = 'GR';
    partType{2} = 'RF';
else
    partLimits = rfLimits;
    partType{1} = 'RF';
end
if rfLimits(2) < pulseSequence.numSteps
    numParts = numParts + 1;
    partType{numParts} = 'RO';
    partLimits = [ partLimits; rfLimits(2)+1, pulseSequence.numSteps ];
end
pulseSequence.numParts   = numParts; % number of parts
pulseSequence.partType   = partType; % type of part
pulseSequence.partLimits = partLimits; % index start/end of part
