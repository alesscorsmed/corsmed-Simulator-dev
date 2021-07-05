function [pulseSequence] = spinEcho2D(...
    acquisition, encoding, mrSystem, expControl)
%
% SEQUENCE.FAMILYSE.SPINECHO2D
%
%	Generates a Standard 2D Spin Echo sequence.
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
functionName = 'sequence.familySE.spinEcho2D';

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
encPlan     = encoding.plan;

%% generate balanced Spin Echo repetition encoding
if acqData.TE > 10e-3 % large TE, use balanced
    [repetition] = sequence.repetitions.balancedFullSpinEcho(...
        acqData, mainRF, refRF, mrSystem, expControl);
else % use unbalanced: allows for shorter TE
    [repetition] = sequence.repetitions.unbalancedFullSpinEcho(...
        acqData, mainRF, refRF, mrSystem, expControl);
end

%% apply spoiling as SWC at the end
if acqData.spoiled == 1
    repetition.swcSignal(end) = 1;
end

%% Apply repetition compression
if expControl.sequence.timeCompression
    repetition = sequence.tools.timeCompression(...
        repetition, expControl);
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
% use the encPlan
peEncodingLevels(:) = encPlan.peEncodings(:);

%% apply encodings
[pulseSequence] = data.simulation.initializeSequence();

% gradients
pulseSequence.gxSignal  = kron(feEncodingLevels,repetition.feSignal);
pulseSequence.gySignal  = kron(peEncodingLevels,repetition.peSignal);
pulseSequence.gzSignal  = kron(seEncodingLevels,repetition.seSignal);
% rf
pulseSequence.rfmSignal = kron(rfEncodingLevels,repetition.rfmSignal);
pulseSequence.rfpSignal = kron(unitRepetition,repetition.rfpSignal);
pulseSequence.rffSignal = kron(unitRepetition,repetition.rffSignal);
% slice selection
pulseSequence.gzSignal  = pulseSequence.gzSignal ...
    + kron(ssEncodingLevels,repetition.ssSignal);
% readouts and swc
pulseSequence.swcSignal = kron(unitRepetition,repetition.swcSignal);
pulseSequence.rxSignal  = kron(rxEncodingLevels,repetition.rxSignal);
pulseSequence.rxSignal(pulseSequence.rxSignal > 0) = 1:numREP*numFE;

%% time
% transform in 2D array and add to each 2nd dimension entry the TR
timeShift = (0:numREP-1)*repetition.time(end);
pulseSequence.time  = repmat(repetition.time,[1,numREP]) + timeShift;
pulseSequence.time  = reshape(pulseSequence.time,[],1);

%% pulse sequences info
pulseSequence.name    	= 'Spin Echo 2D';
pulseSequence.type   	= 'SE2D';
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.totalTime = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps  = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs    = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps   = numREP; % number of total repetitions
pulseSequence.numEnc    = numREP; % number of encoding repetitions
pulseSequence.numShots  = numREP; % number of shots
pulseSequence.numTL     = 1; % number of reps in shot
pulseSequence.TE        = repetition.TE; % repetition TE
pulseSequence.TR        = repetition.TR; % repetition TR
pulseSequence.TI        = 0; % preparation TI
pulseSequence.ESP       = 0; % echo spacing
pulseSequence.effTE     = repetition.TE; % effective echo time
pulseSequence.effTR     = repetition.TR; % effective TR
pulseSequence.effTI     = 0; % effective TI (to center of K space)
% maximum interleaving possible (slices in a TR)
pulseSequence.maxIL     = floor(repetition.TR/repetition.minTR);

%% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts = repetition.numParts*numREP;
% part types
pulseSequence.partType = repmat(repetition.partType.',[numREP,1]).';
% this is a bit involved, but avoids loops...
shift = (0:numREP-1)*repetition.numSteps;
partLimits = repetition.partLimits.';
partLimits = repmat(partLimits(:),[1,numREP]) + shift;
pulseSequence.partLimits = reshape(partLimits,2,[]).';

%% same for the index Limits of the readouts
rxLimitsCol            = repetition.numSteps*(0:numREP-1).';
pulseSequence.rxLimits = [rxLimitsCol + repetition.rxLimits(1), ...
    rxLimitsCol + repetition.rxLimits(2)]; 

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
    fprintf(fid, '\n  TE                  %.3fms', 1e3*pulseSequence.TE);
    fprintf(fid, '\n  TR                  %.3fms', 1e3*pulseSequence.TR);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
    
end

