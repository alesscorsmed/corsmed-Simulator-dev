function [pulseSequence] = sliceHopGRE(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYGRE.SLICEHOPGRE
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
functionName = 'sequence:familySE:sliceHopGRE';

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
effTE = acqData.TE;

if prepIR.Apply
    %% generate IR repetition
    % find TI time to the effTE
    % prepIR.TI = prepIR.TI - effTE;
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

% verify concatenations feasibility
maxConcatenations = floor(acqData.TR/repetitionGRE.TR);
if acqData.concatenations > maxConcatenations
    minTR = ceil(repetitionGRE.TR*1e3*acqData.concatenations);
    msg = sprintf( ['Number of concatenations too large for current TR. ',...
        'Maximum concatenations for current TR is %d. ',...
        'Minimum TR for %d concatenations is %.1fms'],...
        maxConcatenations, acqData.concatenations, minTR );
    if isfield(expControl, 'application') && strcmpi(expControl.application, 'edutool')
        eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    else
        ME = MException(['userwarning:',functionName], '%s', msg);
        throw(ME);
    end
end

%% define number of encodings and repetitions, based on sequence params
numFE   = encPlan.numFE;
numPE   = encPlan.numPE;
numSL   = acqData.numSlices; % number of slices
numCN   = acqData.concatenations; % number of concatenations
% new number of slices to deal with concatenations, some may be zero reps
numSL   = numCN*ceil(numSL/numCN); % multiple of concatenations
% total number of repetitions (including dummy ones)
numREP  = numSL*numPE; 

%% modify the PE and RF Freq or RF Phase
peEncodings = encPlan.peEncodings;
phEncodings = zeros(numSL,1);
ffEncodings = zeros(numSL,1);

%% fill slice encodings
numSlices      = acqData.numSlices;
% sliceThickness = acqData.sliceThickness;
% sliceGap       = acqData.sliceGap;
% slabThickness  = numSlices*sliceThickness + (numSlices-1)*sliceGap;
% % compute centers of slices (decreasing order)
% zPositions     = (sliceThickness-slabThickness)/2 ... % start from bottom
%         + linspace(numSlices-1,0,numSlices)*(sliceThickness + sliceGap);
zPositions     = encPlan.zPositions;
% based on slice selection gradient level, frequency shift to those z
%  minus sign due to use in simulation
zFreqShift     = -acqData.gamma*repetitionGRE.ssGradLevel*zPositions;
% use either phase or frequency
if acqData.usePhaseSS
    % use phase for slice selection: rfPhase = 2*pi*freq*t
    phEncodings(1:numSlices) = 2*pi*zFreqShift(:); % time will be added in rep
else
    % use frequency shift for slice selection: rfFreq
    ffEncodings(1:numSlices) = zFreqShift(:);
end

%% prepare the encoding levels to apply to each rep
feEncodingLevels    = zeros(numSL,numPE); % frequency encoding scaling
peEncodingLevels    = zeros(numSL,numPE); % phase encoding
seEncodingLevels    = zeros(numSL,numPE); % for 3D, the slice phase enc
rxEncodingLevels    = zeros(numSL,numPE); % open or not the RX
% RF levels
rfEncodingLevels    = zeros(numSL,numPE); % to scale the RF
phEncodingLevels    = zeros(numSL,numPE); % to add to the RF Phase
ffEncodingLevels    = zeros(numSL,numPE); % RF frequency shift
ssEncodingLevels    = zeros(numSL,numPE); % slice selection grad

%% fill the changing encoding levels
% note that the kron product (reversed order in PE) is done so that
% for each Phase encoding, all the 3D encodings are applied
peEncodingLevels(:) = kron(peEncodings,ones(numSL,1));
phEncodingLevels(:) = kron(phEncodings,ones(1,numPE));
ffEncodingLevels(:) = kron(ffEncodings,ones(1,numPE));

%% fill rest of data and zero dummy reps (to fill the concatenations)
feEncodingLevels(1:numSlices,:) = 1;
rxEncodingLevels(1:numSlices,:) = 1;
rfEncodingLevels(1:numSlices,:) = 1;
if ~expControl.sequence.deactivateSS
    ssEncodingLevels(1:numSlices,:) = 1;
end

%% prepare repetition combination data
numTL                   = numCN; % each shot with number of concatenations
numSHOT                 = numREP/numTL; % as many shots as required
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
    pulseSequence.name	= 'IR Muli-Slice GRE';
    pulseSequence.type	= 'IR-SH-GRE';
else
    pulseSequence.name	= 'Muli-Slice GRE';
    pulseSequence.type	= 'SH-GRE';
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
pulseSequence.effTR     = acqData.shotTR; % effective TR
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


