function [pulseSequence] = MOLLI(...
    acquisition, encoding, mrSystem, expControl, anatomicalModel)
%
% SEQUENCE.FAMILYSSFP.MOLLI
%
%	Generates a MOLLI sequence.
%
% INPUT
%   acquisition         
%   mrSystem        
%   expControl      
%
% OUTPUT
%   pulseSequence   pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.familySSFP.MOLLI';

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

%% get required info
acqData     = acquisition.data;
mainRF      = acquisition.mainRF;
prepIR      = acquisition.prepIR;
encPlan     = encoding.plan;

%% extract info for MOLLI
cardiacCycleDur_msec = (60/anatomicalModel.HR)*1000;

[schemeMOLLI,pauseMOLLI,TIs] = ...
    sequence.tools.extractMOLLIspecs(acquisition.data.MOLLI,...
    cardiacCycleDur_msec);

%% bSSFP readout
% generate balanced GRE repetition encoding
[repetition] = sequence.repetitions.balancedGradEcho(...
    acqData, mainRF, mrSystem, expControl);

% Apply repetition compression
if expControl.sequence.timeCompression
    repetition = sequence.tools.timeCompression(...
        repetition, expControl);
end

% Apply partial echo if needed
if encPlan.startFE > 1
    repetition = sequence.tools.applyPartialEcho(...
        repetition, encPlan.startFE, encPlan.numCE, expControl);
end

% define number of encodings and repetitions, based on sequence
numFE   = encPlan.numFE;
numPE   = encPlan.numPE;
numSE   = encPlan.numSE;
numENC  = numSE*numPE;

% prepare the encoding levels to apply
numREP              = numENC + acqData.dummySSFP; % add dummy to real enc
feEncodingLevels    = ones(numREP,1); % frequency encoding scaling
peEncodingLevels    = zeros(numREP,1); % phase encoding
seEncodingLevels    = zeros(numREP,1); % for 3D, the slice phase enc
rfEncodingLevels    = ones(numREP,1); % to scale the RF
% to define the slice selection
if expControl.sequence.deactivateSS
    ssEncodingLevels = zeros(numREP,1); % zero slice selection grad
else
    ssEncodingLevels = ones(numREP,1);
end
rxEncodingLevels    = zeros(numREP,1);
unitRepetition      = ones(numREP,1);

% 2D, so we only need to modify the PE
% set PE, RX active only for not dummy reps
idxDummy                        = 1:acqData.dummySSFP;
idxRealReps                     = acqData.dummySSFP+1:numREP;
peEncodingLevels(idxRealReps)   = encPlan.peEncodings(:);
rxEncodingLevels(idxRealReps)   = 1;

% effective TE
% add time of repetitions until center of K-space. Found by minmum encoding 
% level of first shot
[~, zeroIndex] = min(abs(peEncodingLevels(:)));
effTE = repetition.time(end)*(zeroIndex(1)-1);
% add echo time
effTE = effTE + repetition.TE;

% RF scaling for balanced SSFP
switch lower(acqData.prepBSSFP)
    case 'half'
        % first encoding is 1/2 FA, rest FA
        rfEncodingLevels(1) = 0.5;
    otherwise %case 'ramp'
        % ramp-up FA from 10% to 100% of FA during dummy
        rfEncodingLevels(idxDummy) = linspace(0.1,1,acqData.dummySSFP);       
end

% apply phase alternation
% apply phase to real reps
rfEncodingLevels(idxRealReps) = rfEncodingLevels(idxRealReps).* ...
    cos(encPlan.encPhase(:));
% apply corresponding phase of dummy reps
rfEncodingLevels(idxDummy) = rfEncodingLevels(idxDummy).* ...
    cos(encPlan.encPhase(acqData.dummySSFP+1:-1:2));

%% IR pulses
% Adjust TIs till the start of the bSSFP readout
TIs = TIs - effTE;

assert(any(TIs>0),...
    sprintf(['\nERROR: %s : The selected TIs do not fit with the design ',...
    'of the bSSFP readout. \n'],functionName)); % @@@

% generate IR repetition
for iIR = size(TIs,2):-1:1
    prepIR.TI           = TIs(1,iIR);
    prepIR.postRFswc    = 1;
    repetitionIR(iIR)   = sequence.repetitions.prepInvRecovery(...
        acqData, prepIR, mrSystem, expControl);
end

% Apply repetition compression
if expControl.sequence.timeCompression
    repetitionIR = sequence.tools.timeCompression(...
        repetitionIR, expControl);
end

%% allocate space
[pulseSequence] = data.simulation.initializeSequence();
% number of steps
numSteps = sum([repetitionIR(:).numSteps]) + ...
    sum(schemeMOLLI)*numREP*repetition.numSteps + size(schemeMOLLI,2);
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.gxSignal  = zeros(numSteps,1);
pulseSequence.gySignal  = zeros(numSteps,1);
pulseSequence.gzSignal  = zeros(numSteps,1);
pulseSequence.rfmSignal = zeros(numSteps,1);
pulseSequence.rfpSignal = zeros(numSteps,1);
pulseSequence.rffSignal = zeros(numSteps,1);
pulseSequence.swcSignal = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1);

% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts = sum([repetitionIR(:).numParts]) + ...
    sum(schemeMOLLI)*numREP*repetition.numParts;
pulseSequence.partType{pulseSequence.numParts} = []; % allocate space
pulseSequence.partLimits = zeros(pulseSequence.numParts,2);

idx                 = 0;
shot                = 0;
currentTimepoint    = 0;
currentPart         = 0;
for scheme = 1:size(schemeMOLLI,2)    
    for LL = 1:schemeMOLLI(1,scheme)
        
        shot = shot + 1;

        %% indexes for IR and for repetitions
        if LL == 1
            idxIR   = idx + (1:repetitionIR(scheme).numSteps);
            idxRep  = idx + repetitionIR(scheme).numSteps + ...
                (1:numREP*repetition.numSteps);          
        else
            idxRep  = idx + (1:numREP*repetition.numSteps);
        end

        %% apply signals to IR
        if LL == 1
            % gradients
            pulseSequence.gxSignal(idxIR)   = repetitionIR(scheme).feSignal;
            pulseSequence.gySignal(idxIR)   = repetitionIR(scheme).peSignal;
            pulseSequence.gzSignal(idxIR)   = repetitionIR(scheme).seSignal;
            % rf
            pulseSequence.rfmSignal(idxIR)  = repetitionIR(scheme).rfmSignal;
            pulseSequence.rfpSignal(idxIR)  = repetitionIR(scheme).rfpSignal;
            pulseSequence.rffSignal(idxIR)  = repetitionIR(scheme).rffSignal;
            % slice selection
            pulseSequence.gzSignal(idxIR)   = pulseSequence.gzSignal(idxIR) ...
                + repetitionIR(scheme).ssSignal;
            % time
            pulseSequence.time(idxIR)       = currentTimepoint + repetitionIR(scheme).time; % @@@
            currentTimepoint                = max(pulseSequence.time(idxIR));
            
            % define partType for the IR
            pulseSequence.partType(currentPart+(1:repetitionIR(scheme).numParts)) = ...
                repetitionIR(scheme).partType;
            
            % Update partLimits: repetitionIR
            pulseSequence.partLimits(currentPart+(1:repetitionIR(scheme).numParts),:) = ...
                repetitionIR.partLimits;
            
            % Update index
            currentPart = currentPart + repetitionIR(scheme).numParts; % update currentPart
        end

        %% apply encodings to repetitions
        % gradients
        pulseSequence.gxSignal(idxRep)  = kron(feEncodingLevels,repetition.feSignal);
        pulseSequence.gySignal(idxRep)  = kron(peEncodingLevels,repetition.peSignal);
        pulseSequence.gzSignal(idxRep)  = kron(seEncodingLevels,repetition.seSignal);
        % rf
        pulseSequence.rfmSignal(idxRep) = kron(rfEncodingLevels,repetition.rfmSignal);
        pulseSequence.rfpSignal(idxRep) = kron(unitRepetition,repetition.rfpSignal);
        pulseSequence.rffSignal(idxRep) = kron(unitRepetition,repetition.rffSignal);
        
        % slice selection
        pulseSequence.gzSignal(idxRep)  = pulseSequence.gzSignal(idxRep) ...
            + kron(ssEncodingLevels,repetition.ssSignal);
        
        % readouts and swc
        pulseSequence.swcSignal(idxRep) = kron(unitRepetition,repetition.swcSignal);
        pulseSequence.rxSignal(idxRep)  = kron(rxEncodingLevels,repetition.rxSignal);
        
        % receiver - Find the inds of this shot and update its values
        receiverShotInd = pulseSequence.rxSignal(1:idxRep(end),1)>0;
        receiverShotInd(1:end-size(idx+1:idxRep(end),2)) = 0;
        pulseSequence.rxSignal(receiverShotInd) = ...
            (shot-1)*numENC*numFE + (1:numENC*numFE);

        % define partType for the bSSFP readout
        pulseSequence.partType(currentPart+(1:numREP*repetition.numParts)) = ...
            repmat(repetition.partType.',[numREP,1]).';
        
        % Update partLimits: repetition
        shift = idx + (0:numREP-1)*repetition.numSteps;
        partLimits = repetition.partLimits.';
        partLimits = repmat(partLimits(:),[1,numREP]) + shift;
        partLimits = reshape(partLimits,2,[]).';
        pulseSequence.partLimits(currentPart+repetitionIR(scheme).numParts+(1:numREP*repetition.numParts),:) = ...
            partLimits + repetitionIR(scheme).numSteps;
        
        % Update rxLimits: limits for the readouts of repetition
        rxLimitsCol            = idx + repetition.numSteps*(acqData.dummySSFP:numREP-1).';
        pulseSequence.rxLimits = [rxLimitsCol + repetition.rxLimits(1), ...
            rxLimitsCol + repetition.rxLimits(2)] + repetitionIR(scheme).numSteps;
        
        % Update indeces
        currentPart = currentPart + numREP*repetition.numParts; % update currentPart
        idx         = idxRep(end) + 1; % update index for the wait time till the next repetition

        %% time
        % transform in 2D array and add to each 2nd dimension entry the TR
        timeShift                   = (0:numREP-1)*repetition.time(end);
        repTime                     = repmat(repetition.time,[1,numREP]) + timeShift;
        repTime                     = reshape(repTime,[],1);
        pulseSequence.time(idxRep)  = currentTimepoint + repTime;
        currentTimepoint            = max(pulseSequence.time(idxRep));
        
        %% Calculate and add the waiting time till next bSSFP readout
        if LL ~= schemeMOLLI(1,scheme)
            waitTime            = round(cardiacCycleDur_msec/1000) - ...
                numREP*repetition.time(end);
            currentTimepoint    = currentTimepoint + waitTime;
            pulseSequence.time(max(idxRep)+1) = currentTimepoint;
        end
        
    end
    %% Calculate and add the waiting time till next IR
    if scheme ~= size(schemeMOLLI,2)   
        waitTime            = (pauseMOLLI(scheme)+1)*round(cardiacCycleDur_msec/1000) - ...
            numREP*repetition.time(end) - repetitionIR(scheme+1).time(end);
        currentTimepoint    = currentTimepoint + waitTime;
        pulseSequence.time(max(idxRep)+1) = currentTimepoint;
    end
end

% @@@ a/2 at the end of bSSFP readout?
% @@@ Why do I get 7 more timepoints than the ones initially calculated?
% @@@ RF and receiver phases should alternate for the entire pulse sequence 

%% pulse sequences info
pulseSequence.name          = 'MOLLI';
pulseSequence.type          = 'IR-bSSFP';
pulseSequence.endEvent      = 'none'; % indicates what happens at the end
pulseSequence.totalTime     = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps      = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs        = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps       = numREP+1;
pulseSequence.numEnc        = numENC; % number of encoding repetitions
pulseSequence.numShots      = numREP; % number of shots
pulseSequence.numTL         = 1; % number of reps in shot
pulseSequence.TE            = repetition.TE; % repetition TE;
pulseSequence.TR            = repetition.TR; % repetition TR;
pulseSequence.TI            = repetitionIR.TI; % preparation TI;
pulseSequence.ESP           =  0; % echo spacing
pulseSequence.effTE         = effTE; % effective echo time
pulseSequence.effTR         = repetition.TR; % effective TR
pulseSequence.effTI         = TIs; % effective TI (to center of K space)
% maximum interleaving possible (slices in a TR)
pulseSequence.maxIL         = 1; % BSSFP sequences do not allow interleaving

%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Number steps        %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Number readouts     %d', pulseSequence.numRxs);
    fprintf(fid, '\n  Number 2D encodings %d', numPE);
    fprintf(fid, '\n  Number repetitions  %d', pulseSequence.numReps);
    fprintf(fid, '\n  TI                  %.3fms', 1e3*pulseSequence.TI);
    fprintf(fid, '\n  TE                  %.3fms', 1e3*pulseSequence.TE);
    fprintf(fid, '\n  TR                  %.3fms', 1e3*pulseSequence.TR);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
end