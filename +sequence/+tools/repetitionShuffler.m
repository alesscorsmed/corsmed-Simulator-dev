function [pulseSequence] = repetitionShuffler(...
    mainRep, repCombination, prepRep, irRep, expControl)
%
% SEQUENCE.TOOLS.REPETITIONSHUFFLER
%
%	Combines the repetitions into a new repetition.
%   The combination is as follows:
%       irRep + prepRep + comb(mainRep,encMatrix(:,1)) + TRwait
%       irRep + prepRep + comb(mainRep,encMatrix(:,2)) + TRwait
%       ...
%       irRep + prepRep + comb(mainRep,encMatrix(:,end)) + TRwait
%
%   This is useful to perform IR prep multi-shot or ETL sequences
%
% INPUT
%   repetition  original repetition    
%   expControl      
%
% OUTPUT
%   compressed  compressed data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence:tools:repetitionShuffler';

% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFlie,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% collect number of main reps per shot (Train Length, TL) and number of shots
numTL   = repCombination.numTL;
numSHOT = repCombination.numSHOT;
TR      = repCombination.TR;

%% total number of steps
% steps and parts per shot
minTR    = 0.0;
numSteps = 0;
numParts = 0;
numReps  = 0;
if ~isempty(irRep)
    baseIdxIR   = 1:irRep.numSteps;
    numSteps    = irRep.numSteps;
    partIdxIR   = 1:irRep.numParts;
    numParts    = irRep.numParts;
    numReps     = 1;
    minTR       = irRep.totalTime;
end
if ~isempty(prepRep)
    baseIdxPrep = 1:prepRep.numSteps;
    numSteps    = numSteps + prepRep.numSteps;
    partIdxPrep = 1:prepRep.numParts;
    numParts    = numParts + prepRep.numParts;
    numReps     = numReps + 1;
    minTR       = minTR + prepRep.totalTime;
end
mainRepSteps    = mainRep.numSteps*numTL;
baseIdxRep      = 1:mainRepSteps;
numSteps        = numSteps + mainRepSteps + 1;
partIdxRep      = 1:mainRep.numParts*numTL;
numParts        = numParts + mainRep.numParts*numTL;
numReps         = numReps + numTL;
minTR           = minTR + numTL*mainRep.totalTime;
% multiply by number of Shots
numSteps = numSteps*numSHOT;
numParts = numParts*numSHOT;
numReps  = numReps*numSHOT;

%% check TR timing
TRwait = TR - minTR;
if TRwait < 0
    msg = sprintf( ['The selected Repetition Time (TR=%.2fms) is too short. ',...
        'Minimum Repetition Time for the current configuration is TR=%.2fms'],...
        TR*1e3, minTR*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end

%% prepare the encoding levels to apply to each shot
feEncodingLevels    = ones(numTL,1);  % frequency encoding scaling
peEncodingLevels    = zeros(numTL,1); % phase encoding
seEncodingLevels    = zeros(numTL,1); % for 3D, the slice phase enc
rfEncodingLevels    = ones(numTL,1);  % to scale the RF
phEncodingLevels    = zeros(numTL,1);  % phase to add to the RF
ffEncodingLevels    = zeros(numTL,1);  % frequency shift to add to the RF
ssEncodingLevels    = ones(numTL,1);  % slice selection scaling
rxEncodingLevels    = ones(numTL,1);  % eith Rx or not
unitRepetition      = ones(numTL,1);

%% allocate space
[pulseSequence] = data.simulation.initializeSequence();
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.gxSignal  = zeros(numSteps,1);
pulseSequence.gySignal  = zeros(numSteps,1);
pulseSequence.gzSignal  = zeros(numSteps,1);
pulseSequence.rfmSignal = zeros(numSteps,1);
pulseSequence.rfpSignal = zeros(numSteps,1);
pulseSequence.rffSignal = zeros(numSteps,1);
pulseSequence.swcSignal = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1);
% parts
pulseSequence.numParts = numParts;
pulseSequence.partType{pulseSequence.numParts} = []; % allocate
pulseSequence.partLimits = zeros(pulseSequence.numParts,2);
% rx Limits
pulseSequence.rxLimits = zeros(numTL*numSHOT,2);
% maximum interleaving possible (slices in a TR)
pulseSequence.maxIL    = floor(TR/minTR);

%% loop on the shots and assign data
% cumulative timeShift
timeShift = 0;
idxShift  = 0;
partShift = 0;
% loop
for shot = 1:numSHOT

    if ~isempty(irRep)
        %% apply signals to IR
        idxIR = idxShift + baseIdxIR;
        % time
        pulseSequence.time(idxIR)      = irRep.time + timeShift;
        % gradients
        pulseSequence.gxSignal(idxIR)  = irRep.feSignal;
        pulseSequence.gySignal(idxIR)  = irRep.peSignal;
        pulseSequence.gzSignal(idxIR)  = irRep.seSignal;
        % rf
        pulseSequence.rfmSignal(idxIR) = irRep.rfmSignal;
        pulseSequence.rfpSignal(idxIR) = irRep.rfpSignal;
        pulseSequence.rffSignal(idxIR) = irRep.rffSignal;
        % slice selection
        if ~expControl.sequence.deactivateSS
            pulseSequence.gzSignal(idxIR)  = pulseSequence.gzSignal(idxIR) ...
                + irRep.ssSignal;
        end
        % parts
        idxPart = partIdxIR + partShift;
        pulseSequence.partType(idxPart)     = irRep.partType;
        pulseSequence.partLimits(idxPart,:) = irRep.partLimits + idxShift;
        % increase shifts
        timeShift = timeShift + irRep.time(end);
        idxShift  = idxShift  + irRep.numSteps;
        partShift = partShift + irRep.numParts; 
    end
    
    if ~isempty(prepRep)
        %% apply signals to prep
        idxPrep = idxShift + baseIdxPrep;
        % time
        pulseSequence.time(idxPrep)      = prepRep.time + timeShift;
        % gradients
        pulseSequence.gxSignal(idxPrep)  = prepRep.feSignal;
        pulseSequence.gySignal(idxPrep)  = prepRep.peSignal;
        pulseSequence.gzSignal(idxPrep)  = prepRep.seSignal;
        % rf
        pulseSequence.rfmSignal(idxPrep) = prepRep.rfmSignal;
        pulseSequence.rfpSignal(idxPrep) = prepRep.rfpSignal;
        pulseSequence.rffSignal(idxPrep) = prepRep.rffSignal;
        % slice selection
        if ~expControl.sequence.deactivateSS
            pulseSequence.gzSignal(idxPrep)  = pulseSequence.gzSignal(idxPrep) ...
                + prepRep.ssSignal;
        end
        % parts
        idxPart = partIdxPrep + partShift;
        pulseSequence.partType(idxPart)     = prepRep.partType;
        pulseSequence.partLimits(idxPart,:) = prepRep.partLimits + idxShift;
        % increase shifts
        timeShift = timeShift + prepRep.time(end);
        idxShift  = idxShift  + prepRep.numSteps;
        partShift = partShift + prepRep.numParts; 
    end
    
    %% gather encodings for the shot
    rxEncodingLevels(:)  = repCombination.rx(:,shot);
    feEncodingLevels(:)  = repCombination.fe(:,shot);
    peEncodingLevels(:)  = repCombination.pe(:,shot);
    seEncodingLevels(:)  = repCombination.se(:,shot);
    rfEncodingLevels(:)  = repCombination.rf(:,shot);
    phEncodingLevels(:)  = repCombination.ph(:,shot);
    ffEncodingLevels(:)  = repCombination.ff(:,shot);
    ssEncodingLevels(:)  = repCombination.ss(:,shot);
    
    %% main rep: apply encodings to TL repetitions
    idxRep = idxShift + baseIdxRep;
    %% time
    % transform in 2D array 
    % and add to each 2nd dimension entry the rep Time
    repTimeShift = (0:numTL-1)*mainRep.time(end);
    repTime  = repmat(mainRep.time,[1,numTL]) + repTimeShift;
    repTime  = reshape(repTime,[],1);
    pulseSequence.time(idxRep) = repTime + timeShift; % add global shift
    %% gradients
    pulseSequence.gxSignal(idxRep)  = kron(feEncodingLevels,mainRep.feSignal);
    pulseSequence.gySignal(idxRep)  = kron(peEncodingLevels,mainRep.peSignal);
    pulseSequence.gzSignal(idxRep)  = kron(seEncodingLevels,mainRep.seSignal);
    %% rf
    pulseSequence.rfmSignal(idxRep) = kron(rfEncodingLevels,mainRep.rfmSignal);
    %% phase and slice slection
    % slice selection phase
    phaseTime = mainRep.time; 
    phaseTime(~mainRep.rfEntries) = 0; % times of the RF
    rfpShift = kron(phEncodingLevels,phaseTime); % phase shift for slice selection
    % add phase with shift for slice slection
    pulseSequence.rfpSignal(idxRep) = kron(unitRepetition,mainRep.rfpSignal)...
        + rem(rfpShift,2*pi);
    %% frequency slice selection
    pulseSequence.rffSignal(idxRep) = kron(unitRepetition,mainRep.rffSignal) ...
        + kron(ffEncodingLevels,mainRep.rfEntries); % add freq
    %% slice selection
    if ~expControl.sequence.deactivateSS
        pulseSequence.gzSignal(idxRep)  = pulseSequence.gzSignal(idxRep) ...
            + kron(ssEncodingLevels,mainRep.ssSignal);
    end
    % readouts and swc
    pulseSequence.swcSignal(idxRep) = kron(unitRepetition,mainRep.swcSignal);
    pulseSequence.rxSignal(idxRep)  = kron(rxEncodingLevels,mainRep.rxSignal);

    %% for splitting into parts (RF vs Non-RF)
    idxPart = partIdxRep + partShift;
    pulseSequence.partType(idxPart) = repmat(mainRep.partType.',[numTL,1]).';
    % this is a bit involved, but avoids loops...
    % for the repetitions
    repPartshift = (0:numTL-1)*mainRep.numSteps;
    partLimits = mainRep.partLimits.';
    partLimits = repmat(partLimits(:),[1,numTL]) + repPartshift;
    partLimits = reshape(partLimits,2,[]).';
    pulseSequence.partLimits(idxPart,:) = partLimits + idxShift;
    
    %% same for the index Limits of the readouts
    idxRx = (shot-1)*numTL+1:shot*numTL;
    rxLimitsCol = mainRep.numSteps*(0:numTL-1).';
    pulseSequence.rxLimits(idxRx,1) = rxLimitsCol + mainRep.rxLimits(1) + idxShift;
    pulseSequence.rxLimits(idxRx,2) = rxLimitsCol + mainRep.rxLimits(2) + idxShift;
    
    % increase shifts
    partShift = partShift + mainRep.numParts*numTL; 
    idxShift  = idxShift  + mainRep.numSteps*numTL + 1;
    % end time to be compliant with TR
    pulseSequence.time(idxShift) = shot*TR;
    timeShift = shot*TR;
    % last point of the part is the last point of the shot
    pulseSequence.partLimits(idxPart(end),2) = idxShift;

end

%% readouts
pulseSequence.totalTime    = pulseSequence.time(end);
pulseSequence.numSteps     = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs       = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps      = numReps;
pulseSequence.numShots     = numSHOT;
pulseSequence.numTL        = numTL;
pulseSequence.rxSignal(pulseSequence.rxSignal > 0) = 1:pulseSequence.numRxs;
    
%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Train Length        %d', numTL);
    fprintf(fid, '\n  Number of Shots     %d', numSHOT);
    fprintf(fid, '\n  Total # reps        %d', pulseSequence.numReps);
    fprintf(fid, '\n  Total # steps       %d', pulseSequence.numSteps);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

