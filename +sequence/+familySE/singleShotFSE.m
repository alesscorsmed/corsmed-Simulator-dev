function [pulseSequence] = singleShotFSE(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYSE.SINGLESHOTFSE
%
%	Generates a (IR-prep) Single Shot Fast Spin Echo sequence.
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
prepIR      = acquisition.prepIR;
encPlan     = encoding.plan;

%% define number of encodings and repetitions, based on sequence
numPE   = encPlan.numPE;
numSE   = encPlan.numSE;
numENC  = numSE*numPE;
numTL   = numENC;
numSHOT = 1;
numREP  = numTL*numSHOT;

%% generate refocusing (180+Encoding) Spin Echo repetition encoding
[repetition180] = sequence.repetitions.refocusSpinEcho(...
    acqData, refRF, mrSystem, expControl);

% if user forced min TE
if acqData.forceMinTE
    userTE = repetition180.TE;
else
    % user provided TE
    userTE  = acqData.TE;
end

%% generate the preparation block (90+FERW)
% update TE to be that provided by refocusing rep
acqData.TE = repetition180.TE;
% pre-encoding gradient
areaFE = sum(repetition180.feSignal.*...
    [repetition180.time(1); diff(repetition180.time)]);
% generate the RF 90 Prep repetition 
[repetition90] = sequence.repetitions.prepSpinEcho(areaFE,...
    acqData, mainRF, mrSystem, expControl);

%% verify timings
minTE = max(repetition180.TE, repetition90.TE);
if userTE < minTE
    msg = sprintf( ['The selected Echo Time (TE=%.2fms) is too short. ',...
        'Minimum Echo Time for the current configuration is TE=%.2fms'],...
        userTE*1e3, minTE*1e3 );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end

if prepIR.Apply
    %% generate IR repetition
    [repetitionIR] = sequence.repetitions.prepInvRecovery(...
        acqData, prepIR, mrSystem, expControl);
else
    %% no IR: empty repetition and no TI
    repetitionIR = [];
    prepIR.TI = 0;
end

%% Apply repetition compression
if expControl.sequence.timeCompression
    if ~isempty(repetitionIR)
        repetitionIR = sequence.tools.timeCompression(...
            repetitionIR, expControl);
    end
    repetition90 = sequence.tools.timeCompression(...
        repetition90, expControl);
    repetition180 = sequence.tools.timeCompression(...
        repetition180, expControl);
end

%% Apply partial echo if needed
if encPlan.startFE > 1
    repetition180 = sequence.tools.applyPartialEcho(...
        repetition180, encPlan.startFE, encPlan.numCE, expControl);
end

%% prepare the encoding levels to apply to each rep
feEncodingLevels    = zeros(numREP,1); % frequency encoding scaling
peEncodingLevels    = zeros(numREP,1); % phase encoding
seEncodingLevels    = zeros(numREP,1); % for 3D, the slice phase enc
rxEncodingLevels    = zeros(numREP,1); % open or not the RX
% RF levels
rfEncodingLevels    = zeros(numREP,1); % to scale the RF
phEncodingLevels    = zeros(numREP,1); % to add to the RF Phase
ffEncodingLevels    = zeros(numREP,1); % RF frequency shift
ssEncodingLevels    = zeros(numREP,1); % slice selection grad

%% 2D, so we only modify the PE
peEncodingLevels(1:numENC) = encPlan.peEncodings;

%% fill rest of useful data (rest left zero)
feEncodingLevels(1:numENC) = 1;
rxEncodingLevels(1:numENC) = 1;
rfEncodingLevels(1:numENC) = 1;
if ~expControl.sequence.deactivateSS
    ssEncodingLevels(1:numENC) = 1;
end

%% prepare repetition combination data
repCombination.numTL    = numTL;
repCombination.numSHOT  = numSHOT;
repCombination.TR       = acqData.TR;
repCombination.rx       = reshape(rxEncodingLevels,numTL,numSHOT); % readout selection
repCombination.fe       = reshape(feEncodingLevels,numTL,numSHOT); % Freq. encoding scaling
repCombination.pe       = reshape(peEncodingLevels,numTL,numSHOT); % Phase encoding scaling
repCombination.se       = reshape(seEncodingLevels,numTL,numSHOT); % Slice encoding scaling
repCombination.rf       = reshape(rfEncodingLevels,numTL,numSHOT); % RF magnitude scaling
repCombination.ph       = reshape(phEncodingLevels,numTL,numSHOT); % Phase addition
repCombination.ff       = reshape(ffEncodingLevels,numTL,numSHOT); % Frequency addition
repCombination.ss       = reshape(ssEncodingLevels,numTL,numSHOT); % Slice Selection scaling

%% combine into sequence
[pulseSequence] = sequence.tools.repetitionShuffler(...
    repetition180, repCombination, repetition90, repetitionIR, expControl);

%% effective TE
effTE = repetition90.time(end);
% add time of repetitions until center of K-space
%   found by minmum encoding level of first shot
[~, zeroIndex] = min(abs(repCombination.pe(:,1)));
effTE = effTE + repetition180.time(end)*(zeroIndex(1)-1);
% add echo time
effTE = effTE + acqData.TE/2;
% TI for TSE is from center of IR to start of RF90
if prepIR.Apply
    effTI = prepIR.TI + effTE; % actual effective TI to echo center
else
    effTI = 0.0;
end

%% pulse sequences info
if prepIR.Apply
    pulseSequence.name	= 'IR Single Shot FSE';
    pulseSequence.type 	= 'IR-SSFSE';
else
    pulseSequence.name 	= 'Single Shot FSE';
    pulseSequence.type 	= 'SSFSE';
end
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.totalTime = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps  = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs    = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps   = numREP; % number of total repetitions
pulseSequence.numEnc    = numENC; % number of encoding repetitions
pulseSequence.numShots  = numSHOT; % number of shots
pulseSequence.numTL     = numTL; % number of reps in shot
pulseSequence.TE        = repetition180.TE; % repetition TE
pulseSequence.TR        = acqData.TR; % repetition TR
pulseSequence.TI        = prepIR.TI; % preparation TI
pulseSequence.ESP       = 2*repetition180.TE; % echo spacing
pulseSequence.effTE     = effTE; % effective echo time
pulseSequence.effTR     = pulseSequence.totalTime; % effective TR
pulseSequence.effTI     = effTI; % effective TI (to center of K space)

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s sequence, elapsed time %.3fs',...
        functionName, pulseSequence.name, toc(tTotal));
    fprintf(fid, '\n  IRL Time            %.3fs', pulseSequence.totalTime);
    fprintf(fid, '\n  Number steps        %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Number readouts     %d', pulseSequence.numRxs);
    fprintf(fid, '\n  Number encodings    %d', pulseSequence.numEnc);
    fprintf(fid, '\n  Number repetitions  %d', pulseSequence.numReps);
    fprintf(fid, '\n  Effective TI        %.3fms', 1e3*pulseSequence.effTI);
    fprintf(fid, '\n  Effective TE        %.3fms', 1e3*pulseSequence.effTE);
    fprintf(fid, '\n  TI                  %.3fms', 1e3*pulseSequence.TI);
    fprintf(fid, '\n  TE                  %.3fms', 1e3*pulseSequence.TE);
    fprintf(fid, '\n  TR                  %.3fms', 1e3*pulseSequence.TR);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
end
