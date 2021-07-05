function [pulseSequence] = MPRAGE(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYGRE.MPRAGE
%
%	Generates a MP-RAGE sequence.
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
functionName = 'sequence.familySE.MPRAGE';

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

%% define number of encodings and repetitions, based on sequence
numFE   = encPlan.numFE;
numPE   = encPlan.numPE;
numSE   = encPlan.numSE;
numREP  = numSE*numPE;
if strcmpi(encPlan.encPhaseDir, 'sp')
    % will do all 3D encodings for each Phase encoding
    numTL   = numSE;
    numSHOT = numPE;
else
    % will do all Phase encodings for each 3D encoding
    numTL   = numPE;
    numSHOT = numSE;
end

%% generate GRE repetition encoding
[repetitionGRE] = sequence.repetitions.unbalancedGradEcho(...
    acqData, mainRF, mrSystem, expControl);

%% apply spoiling as SWC at the end
if acqData.spoiled == 1
    repetitionGRE.swcSignal(end) = 1;
end

%% no preparation repetition
repetitionPrep = [];

%% effective TE
% time of repetitions until center of K-space in 3D dimension
effTE = repetitionGRE.time(end)*(ceil(numTL/2) -1);
% add echo time
effTE = effTE + acqData.TE/2;

if prepIR.Apply
    %% generate IR repetition
    % find TI time to the effTE
    prepIR.TI = prepIR.TI - effTE;
    effTI = prepIR.TI + effTE; % actual effective TI
    % generate repetition
    [repetitionIR] = sequence.repetitions.prepInvRecovery(...
        acqData, prepIR, mrSystem, expControl);
else
    %% no IR: empty repetition and no TI
    repetitionIR = [];
    prepIR.TI = 0;
    effTI = 0;
end

%% Apply repetition compression
if expControl.sequence.timeCompression
    if ~isempty(repetitionIR)
        repetitionIR = sequence.tools.timeCompression(...
            repetitionIR, expControl);
    end
    repetitionGRE = sequence.tools.timeCompression(...
        repetitionGRE, expControl);
end

%% Apply partial echo if needed
if encPlan.startFE > 1
    repetitionGRE = sequence.tools.applyPartialEcho(...
        repetitionGRE, encPlan.startFE, encPlan.numCE, expControl);
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

%% 3D, so we only modify the PE and SE
peEncodings = encPlan.peEncodings;
seEncodings = encPlan.seEncodings;
if strcmpi(encPlan.encPhaseDir, 'sp')
    % note that the kron product (reversed order in PE) is done so that
    % for each Phase encoding, all the 3D encodings are applied
    peEncodingLevels(:) = kron(peEncodings,ones(numSE,1));
    seEncodingLevels(:) = kron(seEncodings,ones(1,numPE));
else
    % note that the kron product (reversed order in SE) is done so that
    % for each 3D encoding, all the Phase encodings are applied
    peEncodingLevels(:) = kron(peEncodings,ones(1,numSE));
    seEncodingLevels(:) = kron(seEncodings,ones(numPE,1));
end

%% fill rest of useful data (rest left zero)
feEncodingLevels(:,:) = 1;
rxEncodingLevels(:,:) = 1;
rfEncodingLevels(:,:) = 1;
if ~expControl.sequence.deactivateSS
    ssEncodingLevels(:,:) = 1;
end

%% prepare repetition combination data
repCombination.numTL    = numTL;
repCombination.numSHOT  = numSHOT;
repCombination.TR       = acqData.shotTR;
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
    repetitionGRE, repCombination, repetitionPrep, repetitionIR, expControl);

%% pulse sequences info
if prepIR.Apply
    pulseSequence.name	= 'MP-RAGE';
    pulseSequence.type	= 'MP-RAGE';
else
    pulseSequence.name 	= 'Multi-Shot GRE';
    pulseSequence.type 	= 'MS-GRE';
end
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.totalTime = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps  = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs    = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps   = numREP; % number of total repetitions
pulseSequence.numEnc    = numREP; % number of encoding repetitions
pulseSequence.numShots  = numSHOT; % number of shots
pulseSequence.numTL     = numTL; % number of reps in shot
pulseSequence.TE        = repetitionGRE.TE; % repetition TE
pulseSequence.TR        = acqData.TR; % repetition TR
pulseSequence.TI        = prepIR.TI; % preparation TI
pulseSequence.ESP       = 0.0; % echo spacing
pulseSequence.effTE     = effTE; % effective echo time
pulseSequence.effTR     = pulseSequence.totalTime; % effective TR
pulseSequence.effTI     = effTI; % effective TI (to center of K space)
% maximum interleaving possible (slices in a TR)
pulseSequence.maxIL     = 1; % 3D sequences do not allow interleaving

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


