function [kSpace] = assembleKspace( rawData, encodingData, expControl )
%
% RECONSTRUCTION.SIGNAL.ASSEMBLEKSPACE
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
functionName = 'reconstructor.signal.assembleKspace';
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
numFE           = encodingPlan.numFE;
numPE           = encodingPlan.numPE;
numSE           = encodingPlan.numSE;
numCE           = encodingPlan.numCE;
% verify dimensions
assert(numReads==numFE*numCE*numPE*numSE, ...
    sprintf('%s : ERROR -- number of readouts does not match dimensions',...
    functionName) );

%% reshape for transformations
timeSignal  = reshape(timeSignal, numFE, numCE, numPE*numSE, numCoils);

%% assemble the K-space
% allocate full kSpace
kSpace = zeros( encodingPlan.kSizeFE, encodingPlan.kSizePE, ...
    encodingPlan.kSizeSE, numCoils, encodingPlan.kSizeCE );
% local kSpace
partKspace = zeros(numFE, numPE*numSE);
fullKspace = zeros( encodingPlan.kSizeFE,...
    encodingPlan.kSizePE, ...
    encodingPlan.kSizeSE );
% values to process the time signal into K-space
% indexes for partial fourier / parallel
feIdx       = encodingPlan.feIncidence;
peIdx       = encodingPlan.peIncidence;
seIdx       = encodingPlan.seIncidence;
% starting indexes for partial fouirer
feStart     = encodingPlan.startFE-1;
peStart     = encodingPlan.startPE-1;
seStart     = encodingPlan.startSE-1;
% phase and EPI info
rxPhase     = exp( 1j*repmat(encodingPlan.encPhase.', numFE, 1) );
rxReverse   = encodingPlan.rxReverse;

%% loop on contrasts and channels, and assemble correct K-space
for ii = 1:numCE
    for jj = 1:numCoils
        
        %% current k-Space
        partKspace(:,:) = timeSignal(:,ii,:,jj);
        
        % TEMP FIX FOR SBR - encodingPlan.encPhase in sbr is as  long as
        % numPE*numSE*numCE
        if length(encodingPlan.encPhase)==numPE
            partKspace = partKspace.*rxPhase;
        else
            rxPhasePerContrast = rxPhase(((ii-1)*numPE)+1:ii*numPE);
            partKspace = partKspace.*rxPhasePerContrast;
        end
        
        %% apply receiver acquisition phase
        
        
        %% reverse order of readouts if required
        partKspace(:,rxReverse>0) = flipud(partKspace(:,rxReverse>0));
        
        %% shift FE lines by 1 if required
        if encodingPlan.kSpaceShift > 0
            partKspace(2:end,rxReverse==0) = partKspace(1:end-1,rxReverse==0);
            partKspace(1,rxReverse==0) = 0.0;
        end
        
        %% assign to correct positions in full K-space
        if strcmpi(encodingPlan.encPhaseDir, 'sp')
            % Order: for each Phase encoding, all the 3D encodings are applied
            fullKspace(feIdx,peIdx,seIdx) = permute( ...
                reshape(partKspace,numFE,numSE,numPE), [1,3,2]);
        else
            % Order: for each 3D encoding, all the Phase encodings are applied
            fullKspace(feIdx,peIdx,seIdx) = reshape(partKspace,numFE,numPE,numSE);
        end
        
        %% fill the zeros of partial fourier with reverse conjugate
        if encodingPlan.startFE > 1
            if expControl.useOldSequence
                idxSource = 1:feStart;
                idxTarget = encodingPlan.kSizeFE:-1:encodingPlan.kSizeFE-feStart+1;
                fullKspace(idxTarget,end:-1:1,end:-1:1) = conj(fullKspace(idxSource,:,:) );
            else
                idxTarget = 1:feStart;
                idxSource = encodingPlan.kSizeFE:-1:encodingPlan.kSizeFE-feStart+1;
                fullKspace(idxTarget,end:-1:1,end:-1:1) = conj(fullKspace(idxSource,:,:) );
            end
        end
        if encodingPlan.startPE > 1
            if expControl.useOldSequence
                idxSource = 1:peStart;
                idxTarget = encodingPlan.kSizePE:-1:encodingPlan.kSizePE-peStart+1;
                fullKspace(end:-1:1,idxTarget,end:-1:1) = conj(fullKspace(:,idxSource,:) );
            else
                idxTarget = 1:peStart;
                idxSource = encodingPlan.kSizePE:-1:encodingPlan.kSizePE-peStart+1;
                fullKspace(end:-1:1,idxTarget,end:-1:1) = conj(fullKspace(:,idxSource,:) );
            end
        end
        if encodingPlan.startSE > 1
            idxTarget = 1:seStart;
            idxSource = encodingPlan.kSizeSE:-1:encodingPlan.kSizeSE-seStart+1;
            fullKspace(end:-1:1,end:-1:1,idxTarget) = conj(fullKspace(:,:,idxSource) );
        end
        
        %% assign to K-space
        kSpace(:,:,:,jj,ii) = fullKspace(:,:,:);
        
    end
    
    %% for each echo, if alternate flip FE direction
    if encodingPlan.altMultiEcho && mod(ii+1,2)
        kSpace(:,:,:,:,ii) = flipud(kSpace(:,:,:,:,ii));
    end
    
end



%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  # Freq. Enc.   %3d', numFE );
    fprintf(fid, '\n  # Phase Enc.   %3d', numPE);
    fprintf(fid, '\n  # Slice Enc.   %3d', numSE);
    fprintf(fid, '\n  K-Space size   %3d x %d x %d x %d x %d',...
        encodingPlan.kSizeFE, encodingPlan.kSizePE, ...
        encodingPlan.kSizeSE, numCoils, encodingPlan.kSizeCE);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end   
end

