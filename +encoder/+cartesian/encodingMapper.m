function [encodingMap] = encodingMapper( encodingPlan, expControl )
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

%% number of entries in the signal
numFE       = encodingPlan.numFE;
numPE       = encodingPlan.numPE;
numSE       = encodingPlan.numSE;
numReads    = numFE*numPE*numSE;
signalIdx   = 1:numReads; % number of the encoding sample

%% reshape for transformations
signalIdx   = reshape(signalIdx, numFE, numPE*numSE);

% values to process the time signal into K-space
% indexes for partial fourier / parallel
feIdx       = encodingPlan.feIncidence;
peIdx       = encodingPlan.peIncidence;
seIdx       = encodingPlan.seIncidence;
% starting indexes for partial fouirer
feStart     = encodingPlan.startFE-1;
peStart     = encodingPlan.startPE-1;
seStart     = encodingPlan.startSE-1;
% EPI / reverse read info
rxReverse   = encodingPlan.rxReverse;

%% receiver phase
rxPhase     = repmat(encodingPlan.encPhase.', numFE, 1);

%% allocate for the incidence matrices
% partial kSpace
partKspace      = zeros(numFE, numPE*numSE);
% full K space incidence
kSpaceIdx       = zeros( ...
    encodingPlan.xSize, encodingPlan.ySize, encodingPlan.zSize );
% Conjugate K space incidence for partial encodings
kSpaceConjFE    = zeros( ...
    encodingPlan.xSize, encodingPlan.ySize, encodingPlan.zSize );
kSpaceConjPE    = zeros( ...
    encodingPlan.xSize, encodingPlan.ySize, encodingPlan.zSize );
kSpaceConjSE    = zeros( ...
    encodingPlan.xSize, encodingPlan.ySize, encodingPlan.zSize );

%% current k-Space
partKspace(:,:) = signalIdx(:,:);

%% reverse order of readouts if required
partKspace(:,rxReverse>0) = flipud(partKspace(:,rxReverse>0));

%% shift FE lines by 1 if required
if encodingPlan.kSpaceShift > 0
    partKspace(2:end,rxReverse==0) = partKspace(1:end-1,rxReverse==0);
    partKspace(1,rxReverse==0)     = 0.0;
end
    
%% assign to correct positions in full K-space
if strcmpi(encodingPlan.encOrder, 'sp')
    % Order: for each Phase encoding, all the 3D encodings are applied
    kSpaceIdx(feIdx,peIdx,seIdx) = permute( ...
        reshape(partKspace,numFE,numSE,numPE), [1,3,2]);
else
    % Order: for each 3D encoding, all the Phase encodings are applied
    kSpaceIdx(feIdx,peIdx,seIdx) = reshape(partKspace,numFE,numPE,numSE);
end

%% assign to the structure
encodingMap.rxPhase     = reshape(rxPhase,[],1);
encodingMap.kSpaceIdx   = reshape(kSpaceIdx,[],1);

%% fill the zeros of partial fourier with reverse conjugate
if encodingPlan.startFE > 1
    idxTarget = 1:feStart;
    idxSource = encodingPlan.xSize:-1:encodingPlan.xSize-feStart+1;
    kSpaceConjFE(idxTarget,end:-1:1,end:-1:1) = kSpaceIdx(idxSource,:,:);
    encodingMap.conjFE = reshape(kSpaceConjFE,[],1);
else
    encodingMap.conjFE = [];
end
if encodingPlan.startPE > 1
    idxTarget = 1:peStart;
    idxSource = encodingPlan.ySize:-1:encodingPlan.ySize-peStart+1;
    kSpaceConjPE(end:-1:1,idxTarget,end:-1:1) = kSpaceIdx(:,idxSource,:);
    encodingMap.conjPE = reshape(kSpaceConjPE,[],1);
else
    encodingMap.conjPE = [];
end
if encodingPlan.startSE > 1
    idxTarget = 1:seStart;
    idxSource = encodingPlan.zSize:-1:encodingPlan.zSize-seStart+1;
    kSpaceConjSE(end:-1:1,end:-1:1,idxTarget) = kSpaceIdx(:,:,idxSource);
    encodingMap.conjSE = reshape(kSpaceConjSE,[],1);
else
    encodingMap.conjSE = [];
end

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  # Freq. Enc.   %3d', numFE );
    fprintf(fid, '\n  # Phase Enc.   %3d', numPE);
    fprintf(fid, '\n  # Slice Enc.   %3d', numSE);
    fprintf(fid, '\n  K-Space size   %3d x %d x %d',...
        encodingPlan.xSize, encodingPlan.ySize, encodingPlan.zSize);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end   
end

