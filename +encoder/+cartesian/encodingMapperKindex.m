function [encodingMap] = encodingMapperKindex( encodingPlan, expControl )
%
% ENCODER.CARTESIAN.ENCODINGMAPPER
%
%	Generates mapping information to go from time sample to
%   correct K-space assembly.
%
% INPUT
%   encodingPlan   struct with K-space info to assemble the K-space         
%   expControl      
%
% OUTPUT
%   encodingMap   struct with incidences of time samples into K-space
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'encoder.cartesian.encodingMapper';
if (nargin < 2)
    ME = MException('Reconstruction:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% size of the full-blown K-space
kSizeFE     = encodingPlan.kSizeFE;
kSizePE     = encodingPlan.kSizePE;
kSizeSE     = encodingPlan.kSizeSE;
kSizeCE     = encodingPlan.kSizeCE;

%% full K space incidence
kSpaceIdx   = zeros( kSizeFE, kSizePE, kSizeSE, kSizeCE );

%% number of entries in the signal
numFE       = encodingPlan.numFE;
numPE       = encodingPlan.numPE;
numSE       = encodingPlan.numSE;
numCE       = encodingPlan.numCE;
numSamples  = numFE*numPE*numSE*numCE;

%% incidences and other modifications
% values to process the time signal into K-space
% incidences for partial fourier / parallel
feIdx   = encodingPlan.feIncidence;
peIdx 	= encodingPlan.peIncidence;
seIdx 	= encodingPlan.seIncidence;

% starting indexes for partial fouirer
feStart = encodingPlan.startFE-1;
peStart = encodingPlan.startPE-1;
seStart = encodingPlan.startSE-1;

% contrast incidence and number of encodings per contrast
ceIdx               = encodingPlan.ceIncidence;
numEncPerContrast   = encodingPlan.encPerContrast;   

%% number the samples
signalIdx = 1:numSamples; % number of the encoding sample

%% partial kSpace: reshape in samples per readout and encodings
if encodingPlan.multipleSingleshot
    partKspace = reshape(signalIdx,numFE,numPE,numCE,numSE);
    partKspace = permute(partKspace, [1, 2, 4, 3] ); % contrast to the end
else
    
    partKspace = reshape(signalIdx, numFE, numCE, numPE*numSE);

    % reshape to take into account encodings per contrast
    partKspace = reshape(partKspace, numFE, numEncPerContrast, numCE, []);
    partKspace = permute(partKspace, [1, 2, 4, 3] ); % contrast to the end
end

%% reshape depending on the encoding order
if strcmpi(encodingPlan.encPhaseDir, 'sp')
    % Order: for each Phase encoding, all the 3D encodings are applied
    partKspace = permute( ...
        reshape(partKspace,numFE,numSE,numPE,numCE), [1,3,2,4]);
else
    % Order: for each 3D encoding, all the Phase encodings are applied
     partKspace = reshape(partKspace,numFE,numPE,numSE,numCE);
end

%% assign to the full kSpace with incidences
kSpaceIdx(feIdx,peIdx,seIdx,ceIdx) = partKspace(:,:,:,:);

%% receiver phase
if encodingPlan.multipleSingleshot
    rxPhase = repmat(encodingPlan.encPhase.', numFE,1);
    rxPhase = repmat(rxPhase(:),1,numCE);
else
    rxPhase = repmat(encodingPlan.encPhase.', numFE*numCE, 1);
end

%% assign to the structure
encodingMap.rxPhase     = reshape(rxPhase,[],1);
encodingMap.kSpaceIdx   = reshape(kSpaceIdx,[],1);
encodingMap.feStart     = feStart;
encodingMap.peStart     = peStart;
encodingMap.seStart     = seStart;
encodingMap.kSizeFE     = kSizeFE;
encodingMap.kSizePE     = kSizePE;
encodingMap.kSizeSE     = kSizeSE;
encodingMap.kSizeCE     = kSizeCE;

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  # Freq.    Enc.   %3d', numFE );
    fprintf(fid, '\n  # Phase    Enc.   %3d', numPE);
    fprintf(fid, '\n  # Slice    Enc.   %3d', numSE);
    fprintf(fid, '\n  # Contrast Enc.   %3d', numCE);
    fprintf(fid, '\n  K-Space size   %d x %d x %d x %d',...
        kSizeFE, kSizePE, kSizeSE, kSizeCE);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end   
end

