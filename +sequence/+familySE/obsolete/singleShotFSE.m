function [pulseSequence] = singleShotFSE(...
    acquisition, mrSystem, expControl)
%
% SEQUENCE.FAMILYSE.SINGLESHOTFSE
%
%	Generates a Single Shot Fast Spin Echo sequence.
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
functionName = 'sequence.familySE.singleShotFSE';

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
refRF       = acquisition.refRF;
kSpaceInfo  = acquisition.kSpaceInfo;

%% generate refocusing (180+Encoding) Spin Echo repetition encoding
[repetition180] = sequence.repetitions.refocusSpinEcho(...
    acqData, refRF, mrSystem, expControl);

%% generate the preparation block (90+FERW)
% update TE to be that provided by refocusing rep
acqData.TE = repetition180.TE;
% pre-encoding gradient
areaFE = sum(repetition180.feSignal.*...
    [repetition180.time(1); diff(repetition180.time)]);
% generate the RF 90 Prep repetition 
[repetition90] = sequence.repetitions.prepSpinEcho(areaFE,...
    acqData, mainRF, mrSystem, expControl);

%% check TR timing
minTR = repetition90.totalTime + acqData.numPE*repetition180.totalTime;
TRwait = acqData.TR - minTR;
if TRwait < 0
    if ~isempty(expControl.connLocalDB)
        msg = sprintf( ['The selected Repetition Time (TR=%.2fms) is too short. ',...
            'Minimum Repetition Time for the current configuration is TR=%.2fms'],...
            acqData.TR*1e3, minTR*1e3 );
        messagesDB.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error', expControl.experimentID, expControl.pulseqID);
    else
        ME = MException('sequence:wrongRepTime',...
            '%s : %.1fms TR too short (minimum %.1fms)',...
            functionName, acqData.TR*1e3, minTR*1e3);
        throw(ME);
    end
else
    TRwait = max(TRwait,1e-12);
end

%% Apply repetition compression
if expControl.sequence.timeCompression
    repetition90 = sequence.tools.timeCompression(...
        repetition90, expControl);
    repetition180 = sequence.tools.timeCompression(...
        repetition180, expControl);
end

%% define number of encodings and repetitions, based on sequence
numFE   = kSpaceInfo.numFE;
numPE   = kSpaceInfo.numPE;
numSE   = kSpaceInfo.numSE;
numREP  = numSE*numPE;

%% prepare the encoding levels to apply
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
rxEncodingLevels    = ones(numREP,1);
unitRepetition      = ones(numREP,1);

%% 2D, so we only need to modify the PE
% use the kSpaceInfo
peEncodingLevels(:) = kSpaceInfo.peEncodings(:);

%% allocate space
[pulseSequence] = data.pulseSequence.initialize();

numSteps = repetition90.numSteps + numREP*repetition180.numSteps + 1;
pulseSequence.time      = zeros(numSteps,1);
pulseSequence.gxSignal  = zeros(numSteps,1);
pulseSequence.gySignal  = zeros(numSteps,1);
pulseSequence.gzSignal  = zeros(numSteps,1);
pulseSequence.rfmSignal = zeros(numSteps,1);
pulseSequence.rfpSignal = zeros(numSteps,1);
pulseSequence.rffSignal = zeros(numSteps,1);
pulseSequence.swcSignal = zeros(numSteps,1);
pulseSequence.rxSignal  = zeros(numSteps,1);

%% indexes for prep and for repetitions
idxPrep = 1:repetition90.numSteps;
idxRep  = repetition90.numSteps+1:numSteps-1;

%% apply signals to prep
% gradients
pulseSequence.gxSignal(idxPrep)  = repetition90.feSignal;
pulseSequence.gySignal(idxPrep)  = repetition90.peSignal;
pulseSequence.gzSignal(idxPrep)  = repetition90.seSignal;
% rf
pulseSequence.rfmSignal(idxPrep) = repetition90.rfmSignal;
pulseSequence.rfpSignal(idxPrep) = repetition90.rfpSignal;
pulseSequence.rffSignal(idxPrep) = repetition90.rffSignal;
% slice selection
pulseSequence.gzSignal(idxPrep)  = pulseSequence.gzSignal(idxPrep) ...
    + repetition90.ssSignal;
% time
pulseSequence.time(idxPrep)      = repetition90.time;

%% apply encodings to repetitions
% gradients
pulseSequence.gxSignal(idxRep)  = kron(feEncodingLevels,repetition180.feSignal);
pulseSequence.gySignal(idxRep)  = kron(peEncodingLevels,repetition180.peSignal);
pulseSequence.gzSignal(idxRep)  = kron(seEncodingLevels,repetition180.seSignal);
% rf
pulseSequence.rfmSignal(idxRep) = kron(rfEncodingLevels,repetition180.rfmSignal);
pulseSequence.rfpSignal(idxRep) = kron(unitRepetition,repetition180.rfpSignal);
pulseSequence.rffSignal(idxRep) = kron(unitRepetition,repetition180.rffSignal);
% slice selection
pulseSequence.gzSignal(idxRep)  = pulseSequence.gzSignal(idxRep) ...
    + kron(ssEncodingLevels,repetition180.ssSignal);
% readouts and swc
pulseSequence.swcSignal(idxRep) = kron(unitRepetition,repetition180.swcSignal);
pulseSequence.rxSignal(idxRep)  = kron(rxEncodingLevels,repetition180.rxSignal);
pulseSequence.rxSignal(pulseSequence.rxSignal > 0) = 1:numREP*numFE;

%% time
% transform in 2D array and add to each 2nd dimension entry the TR
timeShift = (0:numREP-1)*repetition180.time(end);
repTime  = repmat(repetition180.time,[1,numREP]) + timeShift;
repTime  = reshape(repTime,[],1);
pulseSequence.time(idxRep) = repTime + repetition90.time(end);

% end time to be compliant with TR
pulseSequence.time(numSteps) = pulseSequence.time(numSteps-1) + TRwait;

%% effective TE
effTE = repetition90.time(end);
% add time of repetitions until center of K-space
effTE = effTE + repetition180.time(end)*(numREP - ceil(acqData.numPE/2) -1);
% add echo time
effTE = effTE + acqData.TE;

%% pulse sequences info
pulseSequence.name         = 'Single Shot FSE';
pulseSequence.type         = 'SSFSE';
pulseSequence.endEvent     = 'none'; % indicates what happens at the end
pulseSequence.totalTime    = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps     = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs       = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps      = numREP+1;
pulseSequence.rep          = repetition180;
% Echo time and TR info
pulseSequence.repTE        = acqData.TE;
pulseSequence.effTE        = effTE;
pulseSequence.effTR        = pulseSequence.totalTime;

%% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts = repetition90.numParts + repetition180.numParts*numREP;
% part types
pulseSequence.partType{pulseSequence.numParts} = []; % allocate
pulseSequence.partType(1:repetition90.numParts) = repetition90.partType;
pulseSequence.partType(repetition90.numParts+1:pulseSequence.numParts) = ...
repmat(repetition180.partType.',[numREP,1]).';
% this is a bit involved, but avoids loops...
pulseSequence.partLimits = zeros(pulseSequence.numParts,2);
pulseSequence.partLimits(1:repetition90.numParts,:) = repetition90.partLimits;
% for the repetitions
shift = (0:numREP-1)*repetition180.numSteps;
partLimits = repetition180.partLimits.';
partLimits = repmat(partLimits(:),[1,numREP]) + shift;
partLimits = reshape(partLimits,2,[]).';
pulseSequence.partLimits(repetition90.numParts+1:end,:) = partLimits ...
    + repetition90.numSteps;

%% same for the index Limits of the readouts
rxLimitsCol            = repetition180.numSteps*(0:numREP-1).';
pulseSequence.rxLimits = [rxLimitsCol + repetition180.rxLimits(1), ...
    rxLimitsCol + repetition180.rxLimits(2)] + repetition90.numSteps; 

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Number steps        %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Number readouts     %d', pulseSequence.numRxs);
    fprintf(fid, '\n  Number encodings    %d', numPE);
    fprintf(fid, '\n  Number repetitions  %d', pulseSequence.numReps);
    fprintf(fid, '\n  Repetition TE       %.3fms', 1e3*pulseSequence.repTE);
    fprintf(fid, '\n  Effective  TE       %.3fms', 1e3*pulseSequence.effTE);
    fprintf(fid, '\n  Effective  TR       %.3fms', 1e3*pulseSequence.effTR);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
end


