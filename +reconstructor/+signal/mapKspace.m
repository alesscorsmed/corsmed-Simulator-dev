function [kSpace] = mapKspace( rawData, encodingData, expControl )
%
% RECONSTRUCTION.SIGNAL.MAPKSPACE
%
%	Assembles the K-space from the time signals.
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructor.signal.mapKspace';
if (nargin < 3)
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

%% get the complex time domain signal
timeSignal = rawData.Sx + 1j* rawData.Sy + rawData.noise;

%% get sizes and info
numCoils        = rawData.numCoils;
numReads        = rawData.numReads;
encodingPlan    = encodingData.plan;
encodingMap     = encodingData.map;
numFE           = encodingPlan.numFE;
numPE           = encodingPlan.numPE;
numSE           = encodingPlan.numSE;
numCE           = encodingPlan.numCE;
% verify dimensions
assert(numReads==numFE*numCE*numPE*numSE, ...
    sprintf('%s : ERROR -- number of readouts does not match dimensions',...
    functionName) );

%% reshape for transformations
timeSignal  = reshape(timeSignal, numReads, numCoils);

%% extract mapping info
feStart     = encodingMap.feStart;
peStart     = encodingMap.peStart;
seStart     = encodingMap.seStart;
kSizeFE     = encodingMap.kSizeFE;
kSizePE     = encodingMap.kSizePE;
kSizeSE     = encodingMap.kSizeSE;
kSizeCE     = encodingMap.kSizeCE;
kSpaceIdx   = encodingMap.kSpaceIdx;
rxPhase     = exp(1j*encodingMap.rxPhase);

%% allocate full K-space and reshape incidence
kSpace = zeros(kSizeFE*kSizePE*kSizeSE*kSizeCE, numCoils);
%% loop on channels, and assemble correct K-space
for jj = 1:numCoils
    %% apply receiver phase and store with incidence in full Kspace
    kSpace(kSpaceIdx > 0, jj) = timeSignal(kSpaceIdx(kSpaceIdx>0),jj).*rxPhase;
end

%% reshape to full size
kSpace = reshape(kSpace, kSizeFE, kSizePE, kSizeSE, kSizeCE, numCoils);

%% apply EPI RO reverse correction
if encodingPlan.epiReverse
    %% reverse read direction when needed
    rxReversePE = encodingPlan.peIncidence(encodingPlan.rxReversePE==0);
    kSpace(:,rxReversePE,:,:,:) = flipud(kSpace(:,rxReversePE,:,:,:));
    %% shift FE lines by 1 if required
    if encodingPlan.kSpaceShift > 0
        % circularly shift not reversed columns and zero first entry
        rxReversePE = encodingPlan.peIncidence(encodingPlan.rxReversePE>0);
        kSpace(:,rxReversePE,:,:,:) = circshift(kSpace(:,rxReversePE,:,:,:),1,1);
        kSpace(1,rxReversePE,:,:,:) = 0.0;
    end
end
    
%% fill the zeros of partial fourier with reverse conjugate
if encodingPlan.conjFill
    if feStart > 1
        idxTarget = 1:feStart;
        idxSource = kSizeFE:-1:kSizeFE-feStart+1;
        kSpace(idxTarget,:,:,:,:) = conj(kSpace(idxSource,:,:,:,:));
    end
    if peStart > 1
        idxTarget = 1:peStart;
        idxSource = kSizePE:-1:kSizePE-peStart+1;
        kSpace(:,idxTarget,:,:,:) = conj(kSpace(:,idxSource,:,:,:));
    end
    if seStart > 1
        idxTarget = 1:seStart;
        idxSource = kSizeSE:-1:kSizeSE-seStart+1;
        kSpace(:,:,idxTarget,:,:) = conj(kSpace(:,:,idxSource,:,:));
    end
end
%% alternate the multi-echo
if encodingPlan.altMultiEcho
    kSpace(:,:,:,2:2:end,:) = flipud(kSpace(:,:,:,2:2:end,:));
end



%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  # Freq.    Enc.   %3d', numFE );
    fprintf(fid, '\n  # Phase    Enc.   %3d', numPE);
    fprintf(fid, '\n  # Slice    Enc.   %3d', numSE);
    fprintf(fid, '\n  # Contrast Enc.   %3d', numCE);
    fprintf(fid, '\n  K-Space size   %d x %d x %d x %d x %d',...
        kSizeFE, kSizePE, kSizeSE, kSizeCE, numCoils);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end   
end

