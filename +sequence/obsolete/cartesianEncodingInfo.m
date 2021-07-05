function [kSpaceInfo] = cartesianEncodingInfo(acqData, expControl)
%
% SEQUENCE.TOOLS.CARTESIANENCODINGINFO
%
%	Generates K-space encoding information from the acquisition.
%   This info is used to assemble the correct K-space.
%
% INPUT
%   acquisition         
%   expControl      
%
% OUTPUT
%   kSpaceInfo   struct with K-space info to assemble the K-space
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.tools.cartesianEncodingInfo';

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

%% image size
kSpaceInfo.imSizeX = acqData.matrixX;
kSpaceInfo.imSizeY = acqData.matrixY;
kSpaceInfo.imSizeZ = acqData.matrixZ;

%% encoding factors: size multipliers in the FE and PE directions
kSpaceInfo.factorX = acqData.samplingFactor;
kSpaceInfo.factorY = acqData.foldoverY;

%% k-space size
kSpaceInfo.xSize = kSpaceInfo.factorX*acqData.numFE;
kSpaceInfo.ySize = kSpaceInfo.factorY*acqData.numPE;
if acqData.is3D
    kSpaceInfo.is3D = 1;
    kSpaceInfo.zSize = acqData.numSlices;
else
    kSpaceInfo.is3D = 0;
    kSpaceInfo.zSize = 1;
end
kSpaceInfo.cSize = acqData.numEchoes;

% if sequence is 3D multi-slice (concatenations or slice hopping)
if contains(lower(acqData.pulseSeqFamilyName), 'conc') ... % concatenations
        || contains(lower(acqData.pulseSeqFamilyName), 'sh') % slice hopping
    kSpaceInfo.is3D = 0; % no 3D K-space
end

%% K-space padding factors:
%   if padded size is <= 1, no padding
kSpaceInfo.xPadFactor = acqData.xPadFactor;
kSpaceInfo.yPadFactor = acqData.yPadFactor;
kSpaceInfo.zPadFactor = acqData.zPadFactor;

%% actual number of acquisition samples
% this can be smaller than sizes above
% due to partial Fourier or parallel imaging
kSpaceInfo.numFE = kSpaceInfo.xSize;
kSpaceInfo.numPE = kSpaceInfo.ySize;
% full K-space Phase encodings
peEncodings = linspace(1,-1,kSpaceInfo.numPE);
if acqData.is3D
    kSpaceInfo.numSE = acqData.numSlices;
    % full K-space 3D encodings
    seEncodings = linspace(-1,1,kSpaceInfo.numSE);
else
    kSpaceInfo.numSE = 1;
    seEncodings = 0;
end
kSpaceInfo.numCN = kSpaceInfo.cSize; % number of contrasts: sets of FE in each repetition

%% check parallel imaging or partial fourier acceleration
if ~strcmpi(acqData.parallelImaging, 'no')
    % parallel acceleration, for now only in the Phase encoding
    kSpaceInfo.rFactorPE = acqData.rFactor;
    kSpaceInfo.rFactorSE = 1;
    % no Fourier acceleration
    kSpaceInfo.fFactorFE = 1.0; % partial Fourier, <= 1.0
    kSpaceInfo.fFactorPE = 1.0; % partial Fourier, <= 1.0
    kSpaceInfo.fFactorSE = 1.0; % partial Fourier, <= 1.0
else
    % no parallel acceleration
    kSpaceInfo.rFactorPE = 1;
    kSpaceInfo.rFactorSE = 1;
    % check partial Fourier
    switch lower(acqData.partialFourier)
        case lower('readConjugate')
            % partial Fourier in the Frequency encoding (readout)
            kSpaceInfo.fFactorFE = acqData.fFactor; % partial Fourier, <= 1.0
            kSpaceInfo.fFactorPE = 1.0;
            kSpaceInfo.fFactorSE = 1.0; 
        case lower('phaseConjugate')
            % partial Fourier in the Phase encoding
            kSpaceInfo.fFactorFE = 1.0; 
            kSpaceInfo.fFactorPE = acqData.fFactor; % partial Fourier, <= 1.0
            kSpaceInfo.fFactorSE = 1.0; 
        case lower('sliceConjugate')
            % partial Fourier in the 3D encoding
            kSpaceInfo.fFactorFE = 1.0;
            kSpaceInfo.fFactorPE = 1.0;
            kSpaceInfo.fFactorSE = acqData.fFactor; % partial Fourier, <= 1.0
        otherwise
            kSpaceInfo.fFactorFE = 1.0; 
            kSpaceInfo.fFactorPE = 1.0;
            kSpaceInfo.fFactorSE = 1.0;
    end
end

%% index of starting encoding level depending on Fourier acc.
% frequency encodings
numFE               = round(kSpaceInfo.fFactorFE*kSpaceInfo.numFE);
kSpaceInfo.startFE  = 1 + kSpaceInfo.numFE - numFE;
% phase encodings
numPE               = round(kSpaceInfo.fFactorPE*kSpaceInfo.numPE);
kSpaceInfo.startPE  = 1 + kSpaceInfo.numPE - numPE;
% slice (3D) encodings
numSE               = round(kSpaceInfo.fFactorSE*kSpaceInfo.numSE);
kSpaceInfo.startSE  = 1 + kSpaceInfo.numSE - numSE;

%% incidence arrays: give the position of acquired sample in the k-Space matrix
%   incidence of readouts in Frequency direction of K-space, useful for partial fourier
kSpaceInfo.feIncidence = ...
    reshape(kSpaceInfo.startFE:kSpaceInfo.numFE,[],1); 
% incidence of encoding in Phase direction of K-space: 
% include stride of parallel acc
kSpaceInfo.peIncidence = ...
    reshape(kSpaceInfo.startPE:kSpaceInfo.rFactorPE:kSpaceInfo.numPE,[],1); 
% incidence of encoding in Slice direction of K-space:
% include stride of parallel acc
kSpaceInfo.seIncidence = ...
    reshape(kSpaceInfo.startSE:kSpaceInfo.rFactorSE:kSpaceInfo.numSE,[],1);

%% select the encoding levels using the incidence matrices, from initial
% reduced encodings after acceleration
kSpaceInfo.peEncodings = reshape(peEncodings(kSpaceInfo.peIncidence),[],1);
kSpaceInfo.seEncodings = reshape(seEncodings(kSpaceInfo.seIncidence),[],1);

%% update number of encodings, including partial and parallel acc.
kSpaceInfo.numFE = length(kSpaceInfo.feIncidence);
kSpaceInfo.numPE = length(kSpaceInfo.peIncidence);
kSpaceInfo.numSE = length(kSpaceInfo.seIncidence);

%% acquisition info for reverse RO lines and PE incidence
kSpaceInfo.rxReverse = zeros(kSpaceInfo.numPE*kSpaceInfo.numSE,1);
if contains(lower(acqData.pulseSeqFamilyName), 'epi')
    kSpaceInfo.rxReverse(2:2:end)   = 1; % reverse every other
    kSpaceInfo.peIncidence          = flipud(kSpaceInfo.peIncidence);
end

%% kSpace shift
if strcmpi(acqData.kspaceshift, 'yes')
    kSpaceInfo.kSpaceShift = 1;
else
    kSpaceInfo.kSpaceShift = 0;
end

%% multiple echoes in alternate direction
kSpaceInfo.altMultiEcho = 1;

%% phase of the acquisition for each encoding
kSpaceInfo.rxPhase  = zeros(kSpaceInfo.numPE*kSpaceInfo.numSE,1); 
if contains(lower(acqData.pulseSeqFamilyName), 'bssfp')
    kSpaceInfo.rxPhase(2:2:end) = pi; % change phase for BSSFP
end

%% Normal encoding order: first phase then slice encodings
if contains(lower(acqData.pulseSeqFamilyName), 'rage') ... %MP-RAGE
        || contains(lower(acqData.pulseSeqFamilyName), 'conc') ... % concatenations
        || contains(lower(acqData.pulseSeqFamilyName), 'sh') % slice hopping
    kSpaceInfo.encOrder = 'SP'; % change order, first slice, then phase
else
    kSpaceInfo.encOrder = 'PS'; % first phase, then slice
end

%% reorder encodings for TSE sequences
if contains(lower(acqData.pulseSeqFamilyName), 'tse') && (acqData.ETL > 1)
%     % do not allow non-multiples
%     if mod(numPE,acqData.ETL)
%         if ~isempty(expControl.connLocalDB)
%             msg = sprintf( ['The number of encodings (%d) is not a divisible',...
%                 'by the selected Echo Train Length (ETL=%d)'],...
%                 numPE, acqData.ETL );
%             messagesDB.errorAndDBentry(expControl.connLocalDB, msg, ...
%                 'cancelled-error', expControl.experimentID, expControl.pulseqID);
%         else
%             ME = MException('sequence:wrongETL',...
%                 '%s : Encoding number %d is no multiple of ETL %d',...
%                 functionName, numPE, acqData.ETL);
%             throw(ME);
%         end
%     else
%         % reorder the encodings in groups of non-consecutive ETL entries
%         originIdx = reshape(1:numPE,[],acqData.ETL).';
%         kSpaceInfo.peIncidence(:) = kSpaceInfo.peIncidence(originIdx);
%         kSpaceInfo.peEncodings(:) = kSpaceInfo.peEncodings(originIdx);
%     end
     
    numShot = ceil(numPE/acqData.ETL);
    idxEnc  = zeros(numPE,1);
    counter = 0;
    for shot = 1:numShot
        % get the encoding indexes corresponding to this shot
        idxLocalEnc = shot:numShot:numPE;
        numLocalEnc = length(idxLocalEnc);
        % re-order shot encodings
        switch lower(acqData.encOrder)
            case 'centric'
                % order phase encoding and incidence in ascending order
                [~,shotOrder] = sort( ...
                    abs(kSpaceInfo.peEncodings(idxLocalEnc)), 'ascend');
            case 'edge'
                % order phase encoding and incidence in descencing order
                [~,shotOrder] = sort( ...
                    abs(kSpaceInfo.peEncodings(idxLocalEnc)), 'descend');
            otherwise
                % natural order from negative to positive
                [~,shotOrder] = sort( ... 
                    kSpaceInfo.peEncodings(idxLocalEnc), 'ascend');
        end
        % store ordered encodings
        idxEnc(counter+1:counter+numLocalEnc) = idxLocalEnc(shotOrder);
        counter = counter + numLocalEnc;
    end
    % apply the reorder to the encoding and incidence
    kSpaceInfo.peEncodings(:) = kSpaceInfo.peEncodings(idxEnc);
    kSpaceInfo.peIncidence(:) = kSpaceInfo.peIncidence(idxEnc);

else
    
    switch lower(acqData.encOrder)
        case 'centric'
            % order phase encoding and incidence in ascending order
            [~,idxEnc] = sort(abs(kSpaceInfo.peEncodings),'ascend');
        case 'edge'
            % order phase encoding and incidence in descencing order
            [~,idxEnc] = sort(abs(kSpaceInfo.peEncodings),'descend');
        otherwise
            % natural order from negative to positive
            [~,idxEnc] = sort(kSpaceInfo.peEncodings,'ascend');
    end
    % do not apply to EPI sequences
    if ~contains(lower(acqData.pulseSeqFamilyName), 'epi')
        kSpaceInfo.peEncodings(:) = kSpaceInfo.peEncodings(idxEnc);
        kSpaceInfo.peIncidence(:) = kSpaceInfo.peIncidence(idxEnc);
    end
end

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
