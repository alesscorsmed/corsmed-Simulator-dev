function [pulseSequence] = singleRF(mainRF, mrSystem, expControl)
%
% SEQUENCE.DUMMY.SINGLERF
%
%	Generates a RFBlock with optional slice selection,
%   and pre/post gradient rewinds.
%
% INPUT
%   duration        total RF duration, in s
%   cycles          number of cycles of the sinc
%   angle           Flip angle, in degrees
%   tstep           time discretization
%   gamma           gyromagnetic ratio
%
% OUTPUT
%   pulseSequence        pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.dummy.singleRF';

if (nargin<3) || isempty(mainRF) || isempty(mrSystem) || isempty(expControl)
    ME = MException('sequence:wrongArgument',...
        '%s : invalid arguments',functionName);
    throw(ME);
end

% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% generate the waveforms for the RF
gamma = 42.577478518e6;
[time,rfm,rfp,rff,gr,rfLimits] = sequence.excitations.sliceSelectRF(...
    mainRF, mrSystem.maxGStrenght, mrSystem.SlewRate,...
    gamma, expControl.sequence.dtRF, expControl );

%% sequences data
[pulseSequence] = data.pulseSequence.initialize();
pulseSequence.type         = 'Single RF';
pulseSequence.endEvent     = 'none'; % indicates what happens at the end
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

% for area compressed sequences
pulseSequence.hArea     = zeros(pulseSequence.numSteps,1); % time deltas for areas
pulseSequence.gxArea    = zeros(pulseSequence.numSteps,1); % x gradient area delta
pulseSequence.gyArea    = zeros(pulseSequence.numSteps,1); % y gradient area delta
pulseSequence.gzArea    = zeros(pulseSequence.numSteps,1); % z gradient area delta

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
    partType{numParts} = 'GR';
    partLimits = [ partLimits; rfLimits(2)+1, pulseSequence.numSteps ];
end
pulseSequence.numParts   = numParts; % number of parts
pulseSequence.partType   = partType; % type of part
pulseSequence.partLimits = partLimits; % index start/end of part

% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  IRL Time         %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Number steps     %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Number readouts  %d', pulseSequence.numRxs);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end