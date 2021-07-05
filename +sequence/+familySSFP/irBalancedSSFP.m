function [pulseSequence] = irBalancedSSFP(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYSSFP.IRBALANCEDSSFP
%
%	Generates a IR bSSFP sequence.
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
functionName = 'sequence.familySSFP.irBalancedSSFP';

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

%% generate IR repetition
[repetitionIR] = sequence.repetitions.prepInvRecovery(...
    acqData, prepIR, mrSystem, expControl);

%% generate balanced GRE repetition encoding
[repetition] = sequence.repetitions.balancedGradEcho(...
    acqData, mainRF, mrSystem, expControl);

%% Apply repetition compression
if expControl.sequence.timeCompression
    repetition = sequence.tools.timeCompression(...
        repetition, expControl);
    repetitionIR = sequence.tools.timeCompression(...
        repetitionIR, expControl);
end

%% Apply partial echo if needed
if encPlan.startFE > 1
    repetition = sequence.tools.applyPartialEcho(...
        repetition, encPlan.startFE, encPlan.numCE, expControl);
end

%% define number of encodings and repetitions, based on sequence
numFE   = encPlan.numFE;
numPE   = encPlan.numPE;
numSE   = encPlan.numSE;
numENC  = numSE*numPE;

%% prepare the encoding levels to apply
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

%% 2D, so we only need to modify the PE
% set PE, RX active only for not dummy reps
idxDummy                        = 1:acqData.dummySSFP;
idxRealReps                     = acqData.dummySSFP+1:numREP;
peEncodingLevels(idxRealReps)   = encPlan.peEncodings(:);
rxEncodingLevels(idxRealReps)   = 1;

%% RF scaling for balanced SSFP
switch lower(acqData.prepBSSFP)
    case 'half'
        % first encoding is 1/2 FA, rest FA
        rfEncodingLevels(1) = 0.5;
    otherwise %case 'ramp'
        % ramp-up FA from 10% to 100% of FA during dummy
        rfEncodingLevels(idxDummy) = linspace(0.1,1,acqData.dummySSFP);       
end

%% apply phase alternation
% apply phase to real reps
rfEncodingLevels(idxRealReps) = rfEncodingLevels(idxRealReps).* ...
    cos(encPlan.encPhase(:));
% apply corresponding phase of dummy reps
rfEncodingLevels(idxDummy) = rfEncodingLevels(idxDummy).* ...
    cos(encPlan.encPhase(acqData.dummySSFP+1:-1:2));

%% apply arbitrary RF scaling (sinusoidal train-like)
if isfield(expControl.sequence,'rfCosTrain') && expControl.sequence.rfCosTrain
    rfEncodingLevels(idxRealReps) = rfEncodingLevels(idxRealReps).* ...
        abs(sin(6*(idxRealReps-idxRealReps(1))*pi/180))';
end

%% allocate space
[pulseSequence] = data.simulation.initializeSequence();
% number of steps
numSteps = repetitionIR.numSteps + numREP*repetition.numSteps;
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.gxSignal  = zeros(numSteps,1);
pulseSequence.gySignal  = zeros(numSteps,1);
pulseSequence.gzSignal  = zeros(numSteps,1);
pulseSequence.rfmSignal = zeros(numSteps,1);
pulseSequence.rfpSignal = zeros(numSteps,1);
pulseSequence.rffSignal = zeros(numSteps,1);
pulseSequence.swcSignal = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1);

%% indexes for IR and for repetitions
idxIR   = 1:repetitionIR.numSteps;
idxRep  = repetitionIR.numSteps+1:numSteps;

%% apply signals to prep
% gradients
pulseSequence.gxSignal(idxIR)  = repetitionIR.feSignal;
pulseSequence.gySignal(idxIR)  = repetitionIR.peSignal;
pulseSequence.gzSignal(idxIR)  = repetitionIR.seSignal;
% rf
pulseSequence.rfmSignal(idxIR) = repetitionIR.rfmSignal;
pulseSequence.rfpSignal(idxIR) = repetitionIR.rfpSignal;
pulseSequence.rffSignal(idxIR) = repetitionIR.rffSignal;
% slice selection
pulseSequence.gzSignal(idxIR)  = pulseSequence.gzSignal(idxIR) ...
    + repetitionIR.ssSignal;
% time
pulseSequence.time(idxIR)      = repetitionIR.time;

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
pulseSequence.rxSignal(pulseSequence.rxSignal > 0) = 1:numENC*numFE;

%% time
% transform in 2D array and add to each 2nd dimension entry the TR
timeShift = (0:numREP-1)*repetition.time(end);
repTime  = repmat(repetition.time,[1,numREP]) + timeShift;
repTime  = reshape(repTime,[],1);
pulseSequence.time(idxRep) = repTime + repetitionIR.time(end);

%% effective TE
% add time of repetitions until center of K-space
%   found by minmum encoding level of first shot
[~, zeroIndex] = min(abs(peEncodingLevels(:)));
effTE = repetition.time(end)*(zeroIndex(1)-1);
% add echo time
effTE = effTE + repetition.TE;
% TI for TSE is from center of IR to start of RF90
effTI = repetitionIR.TI + effTE; % actual effective TI to echo center

%% pulse sequences info
pulseSequence.name         = 'IR Balanced SSFP';
pulseSequence.type         = 'IR-bSSFP';
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.totalTime = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps  = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs    = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps   = numREP+1;
pulseSequence.numEnc    = numENC; % number of encoding repetitions
pulseSequence.numShots  = numREP; % number of shots
pulseSequence.numTL     = 1; % number of reps in shot
pulseSequence.TE        = repetition.TE; % repetition TE;
pulseSequence.TR        = repetition.TR; % repetition TR;
pulseSequence.TI        = repetitionIR.TI; % preparation TI;
pulseSequence.ESP       =  0; % echo spacing
pulseSequence.effTE     = effTE; % effective echo time
pulseSequence.effTR     = repetition.TR; % effective TR
pulseSequence.effTI     = effTI; % effective TI (to center of K space)
% maximum interleaving possible (slices in a TR)
pulseSequence.maxIL     = 1; % BSSFP sequences do not allow interleaving

%% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts = repetitionIR.numParts + repetition.numParts*numREP;
% part types
pulseSequence.partType{pulseSequence.numParts} = []; % allocate
pulseSequence.partType(1:repetitionIR.numParts) = repetitionIR.partType;
pulseSequence.partType(repetitionIR.numParts+1:pulseSequence.numParts) = ...
    repmat(repetition.partType.',[numREP,1]).';
% this is a bit involved, but avoids loops...
pulseSequence.partLimits = zeros(pulseSequence.numParts,2);
pulseSequence.partLimits(1:repetitionIR.numParts,:) = repetitionIR.partLimits;
% for the repetitions
shift = (0:numREP-1)*repetition.numSteps;
partLimits = repetition.partLimits.';
partLimits = repmat(partLimits(:),[1,numREP]) + shift;
partLimits = reshape(partLimits,2,[]).';
pulseSequence.partLimits(repetitionIR.numParts+1:end,:) = partLimits ...
    + repetitionIR.numSteps;

%% same for the index Limits of the readouts
rxLimitsCol            = repetition.numSteps*(acqData.dummySSFP:numREP-1).';
pulseSequence.rxLimits = [rxLimitsCol + repetition.rxLimits(1), ...
    rxLimitsCol + repetition.rxLimits(2)] + repetitionIR.numSteps; 

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


