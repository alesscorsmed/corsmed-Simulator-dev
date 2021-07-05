function [pulseSequence] = turboSE(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYSE.TURBOSE
%
%	Generates a (IR-prepared) Turbo Spin Echo sequence.
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
functionName = 'sequence.familySE.turboSE';

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
numFE   = encPlan.numFE;
numPE   = encPlan.numPE;
numSE   = encPlan.numSE;
numENC  = numSE*numPE;
numSHOT = ceil(numENC/acqData.ETL);
numTL   = ceil(numENC/numSHOT);
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
    % generate repetition
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

%% prepare repetition combination data
repCombination.numTL    = numTL;
repCombination.numSHOT  = numSHOT;
repCombination.TR       = acqData.TR;
repCombination.rx       = zeros(numTL,numSHOT); % readout selection
repCombination.fe       = zeros(numTL,numSHOT); % Freq. encoding scaling
repCombination.pe       = zeros(numTL,numSHOT); % Phase encoding scaling
repCombination.se       = zeros(numTL,numSHOT); % Slice encoding scaling
repCombination.rf       = zeros(numTL,numSHOT); % RF magnitude scaling
repCombination.ph       = zeros(numTL,numSHOT); % Phase addition
repCombination.ff       = zeros(numTL,numSHOT); % Frequency addition
repCombination.ss       = zeros(numTL,numSHOT); % Slice Selection scaling

%% Prepare 3D encoding, so we only modify the PE and SE
% note that the kron product (reversed order in SE) is done so that
% for each 3D encoding, all the phase encodings are applied
peEncodings         = encPlan.peEncodings;
seEncodings         = encPlan.seEncodings;
peEncodingLevels    = kron(peEncodings,ones(1,numSE));
seEncodingLevels    = kron(seEncodings,ones(numPE,1));  

%% since there may be shots with smaller TL, loop and assign each Shot
counter = 0;
for shot = 1:numSHOT
    % get the number of reps in the shot
    numLocalEnc = encPlan.encPerShot(shot); % number of encodings in the shot
    % get current encoding indexes
    idxEnc = counter+1:counter+numLocalEnc;
    % assign current encodings
    repCombination.rx(1:numLocalEnc,shot) = 1; % readout selection
    repCombination.fe(1:numLocalEnc,shot) = 1; % Freq. encoding scaling
    repCombination.pe(1:numLocalEnc,shot) = peEncodingLevels(idxEnc); % Phase encoding scaling
    repCombination.se(1:numLocalEnc,shot) = seEncodingLevels(idxEnc); % Slice encoding scaling
    repCombination.rf(1:numLocalEnc,shot) = 1; % RF magnitude scaling
    repCombination.ph(1:numLocalEnc,shot) = 1; % Phase addition
    repCombination.ff(1:numLocalEnc,shot) = 1; % Frequency addition
    repCombination.ss(1:numLocalEnc,shot) = 1; % Slice Selection scaling
    % increase counter
    counter = counter + numLocalEnc;
end
% for slice encoding: active or not
if expControl.sequence.deactivateSS
    repCombination.ss(:) = 0;
end

%% combine into sequence
[pulseSequence] = sequence.tools.repetitionShuffler(...
    repetition180, repCombination, repetition90, repetitionIR, expControl);

%% effective TE
effTE = repetition90.time(end);
% add time of repetitions until center of K-space
%   found by minmum encoding level of first shot
[~, zeroIndex] = min(repCombination.pe( (repCombination.pe(:,1) >=0), 1));
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
    pulseSequence.name	= 'IR Turbo SE';
    pulseSequence.type	= 'IR-TSE';
else
    pulseSequence.name	= 'Turbo SE';
    pulseSequence.type	= 'TSE';
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
% if 3D force maximum interleaving possible (slices in a TR) to 1
if numSE > 1
    pulseSequence.maxIL     = 1; % 3D sequences do not allow interleaving
end

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


