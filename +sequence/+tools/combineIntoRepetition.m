function [repetition] = combineIntoRepetition( ...
    rfTime, rfm, rfp, rff, grSlice, rfLimits, ...
    grTime, signalFE, signalPE, signalSE, signalRX, rxLimits, ...
    midTime, repTime )
%
% SEQUENCE.TOOLS.COMBINEINTOREPETITION
%
%	combines waveforms into a repetition with given timing.
%
% INPUT
%
% OUTPUT
%   repetition  assembled repetition    
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:tools:combineIntoRepetition';

%% correct RW times to generate correct TE
grTime  = grTime + midTime;

%% append midtime as part of the encoding, with zero signals
grTime   = [midTime; grTime(:)];
signalFE = [0.0; signalFE(:)];
signalPE = [0.0; signalPE(:)];
signalSE = [0.0; signalSE(:)];
signalRX = [0.0; signalRX(:)];
rxLimits = rxLimits + 1; % added an extra point at beginning

%% combine into the repetition: unify times
repetition.time = reshape(union(rfTime, grTime),[],1);
% increase the number to match the required TR
numSteps = length(repetition.time);
if repTime > repetition.time(end)
    numSteps = numSteps + 1;
    repetition.time(numSteps,1) = repTime;
end
% time differences
repetition.timeDiff = reshape(diff([0; repetition.time]),[],1);

%% allocate waveforms
repetition.rxSignal  = zeros(numSteps,1); % receiver readout
repetition.swcSignal = zeros(numSteps,1); % software crusher
repetition.feSignal  = zeros(numSteps,1); % freq encoding (x) gradient
repetition.peSignal  = zeros(numSteps,1); % phase encoding (y) gradient
repetition.seSignal  = zeros(numSteps,1); % 3D phase (slice) encoding (z) gradient
repetition.rfmSignal = zeros(numSteps,1); % RF magnitude
repetition.rfpSignal = zeros(numSteps,1); % RF phase
repetition.rffSignal = zeros(numSteps,1); % RF frequency
repetition.ssSignal  = zeros(numSteps,1); % slice selection gradient
repetition.rfEntries = zeros(numSteps,1); % ones or zeros depending on RF

%% find incidence of RF and ENC times into new time array
% get indexes of RF arrays into new time array
[~,~,rfIndex] = intersect(rfTime,repetition.time,'stable');
% get indexes of RF arrays into new time array
[~,~,grIndex] = intersect(grTime,repetition.time,'stable');

%% find intertwinned entries (if any)
commonIndex = grIndex(1):rfIndex(end);
% assign overlapping signals: make sure to keep area
if ~isempty(commonIndex)
    
    %% interpolate the Encoding signals
    % time points before the starting of the encoding:
    iQuery = commonIndex <= grIndex(1);
    %   zero them
    repetition.feSignal(commonIndex(iQuery)) = 0.0;
    repetition.peSignal(commonIndex(iQuery)) = 0.0;
    repetition.seSignal(commonIndex(iQuery)) = 0.0;
    % time points after the starting of the encoding
    iQuery = commonIndex > grIndex(1);
    %   interpolate using next neighbor: keep areas 
    repetition.feSignal(commonIndex(iQuery)) = interp1( grTime, signalFE,...
        repetition.time(commonIndex(iQuery)), 'next');
    repetition.peSignal(commonIndex(iQuery)) = interp1( grTime, signalPE,...
        repetition.time(commonIndex(iQuery)), 'next');
    repetition.seSignal(commonIndex(iQuery)) = interp1( grTime, signalSE,...
        repetition.time(commonIndex(iQuery)), 'next');
    % time points after the end of the encoding:
    iQuery = commonIndex >= grIndex(end);
    %   zero them
    repetition.feSignal(commonIndex(iQuery)) = 0.0;
    repetition.peSignal(commonIndex(iQuery)) = 0.0;
    repetition.seSignal(commonIndex(iQuery)) = 0.0;

    %% interpolate the SS gradient points
    % time points before the end of the ss:
    iQuery = commonIndex <= rfIndex(end);
    %   interpolate using next neighbor: keep areas 
    repetition.ssSignal(commonIndex(iQuery)) = interp1( rfTime, grSlice,...
        repetition.time(commonIndex(iQuery)), 'next');
    % time points after the end of the ss:
    iQuery = commonIndex > rfIndex(end);
    %   zero them
    repetition.ssSignal(commonIndex(iQuery)) = 0.0;

end

%% assign RF signals
% assign signal entries to corrsponging positions
repetition.rfmSignal(rfIndex)   = rfm;
repetition.rfpSignal(rfIndex)   = rfp;
repetition.rffSignal(rfIndex)   = rff;
repetition.ssSignal(rfIndex)    = grSlice;
% correct the rfLimits (if needed) and define the actual RF entries
rfLimits(:) = rfIndex(rfLimits(:));
repetition.rfEntries(rfLimits(1):rfLimits(2)) = 1;

%% assign GR signals
% assign signal entries to corrsponging positions
repetition.feSignal(grIndex) = signalFE(:);
repetition.peSignal(grIndex) = signalPE(:);
repetition.seSignal(grIndex) = signalSE(:);
repetition.rxSignal(grIndex) = signalRX(:);
% correct the position of the rxLimits
rxLimits(:) = grIndex(rxLimits(:));

%% pulse sequences info
repetition.type         = 'General';
repetition.totalTime    = repetition.time(end); % total time in seconds
repetition.numSteps     = length(repetition.time); % number of time steps
repetition.numRxs       = nnz(repetition.rxSignal); % number of readout points
% slice selection gradient Amplitude
repetition.ssGradLevel  = 0;

%% for splitting into parts (RF vs Non-RF)
% multiple parts, depending on its characteristics
numParts = 1;
% RF
if rfLimits(1) > 1
    numParts = numParts + 1;
    partLimits = [1, rfLimits(1)-1; rfLimits];
    partType{1} = 'GR';
    partType{2} = 'RF';
else
    partLimits = rfLimits;
    partType{1} = 'RF';
end
% rest may have readouts
numParts = numParts + 1;
partLimits = [partLimits; rfLimits(2)+1, repetition.numSteps];
if repetition.numRxs > 0
    partType{numParts} = 'RO';
else
    partType{numParts} = 'GR';
end
% assign
repetition.numParts   = numParts; % number of parts
repetition.partType   = partType; % type of part
repetition.partLimits = partLimits; % index start/end of part
repetition.rxLimits   = rxLimits; % index start/end of readouts
% update TE / TR
repetition.TE         = 0;
repetition.TR         = repetition.time(end);



