function [encodingPlan] = encodingPlanner(acqData, expControl)
%
% ENCODER.CARTESIAN.ENCODINGPLANNER
%
%	Generates K-space encoding information from the acquisition.
%   This info is used to assemble the correct K-space.
%
% INPUT
%   acquisition         
%   expControl      
%
% OUTPUT
%   encodingPlan   struct with K-space info to assemble the K-space
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'encoder.cartesian.encodingPlanner';

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

%% image size and FOV
encodingPlan.imSizeX = acqData.matrixX;
encodingPlan.imSizeY = acqData.matrixY;
encodingPlan.imSizeZ = acqData.matrixZ;
encodingPlan.imFovX  = acqData.fovFE;
encodingPlan.imFovY  = acqData.fovPE;
encodingPlan.imFovZ  = acqData.fovSE;

%% encoding factors: size multipliers in the FE and PE directions
encodingPlan.factorX = acqData.samplingFactorFE;
encodingPlan.factorY = acqData.samplingFactorPE;
encodingPlan.factorZ = acqData.samplingFactorSE;

%% k-space size
encodingPlan.kSizeFE = encodingPlan.factorX*acqData.numFE;
encodingPlan.kSizePE = encodingPlan.factorY*acqData.numPE;
if acqData.is3D
    encodingPlan.is3D = 1;
    encodingPlan.kSizeSE = encodingPlan.factorZ*acqData.numSlices;
else
    encodingPlan.is3D = 0;
    encodingPlan.kSizeSE = 1;
end
encodingPlan.kSizeCE = acqData.numEchoes;

% (relative) z coordinates of each slice after rotation / translation
encodingPlan.zPositions = acqData.zSliceCoord;
% check the slice direction to determin if the order is reversed due to rot
encodingPlan.reversedSlices = 0;
if length(encodingPlan.zPositions) > 1 ...
        && encodingPlan.zPositions(1) > encodingPlan.zPositions(2)
    % User played with rotations, slice order reversed (going down)
    % full K-space 3D encodings: descending
    encodingPlan.reversedSlices = 1;
end

% if sequence is 3D multi-slice (concatenations or slice hopping)
if contains(lower(acqData.pulseSeqFamilyName), 'conc') ... % concatenations
        || contains(lower(acqData.pulseSeqFamilyName), 'sh') % slice hopping
    encodingPlan.is3D = 0; % no 3D K-space
    % modify the image Size and FOV in Z to interpole the maps
    encodingPlan.imSizeZ = acqData.numSlices;
    encodingPlan.imFovZ  = acqData.numSlices * ...
        ( acqData.sliceThickness + acqData.sliceGap )...
        - acqData.sliceGap;
end

%% K-space padding factors:
%   if padded size is <= 1, no padding
encodingPlan.xPadFactor = acqData.xPadFactor;
encodingPlan.yPadFactor = acqData.yPadFactor;
encodingPlan.zPadFactor = acqData.zPadFactor;

%% actual number of acquisition samples
% this can be smaller than sizes above
% due to partial Fourier or parallel imaging
encodingPlan.numFE = encodingPlan.kSizeFE;
encodingPlan.numPE = encodingPlan.kSizePE;
% full K-space Phase encodings
peEncodings = linspace(1,-1,encodingPlan.numPE);
if acqData.is3D
    encodingPlan.numSE = acqData.numSlices;
    if encodingPlan.numSE > 1
        if encodingPlan.reversedSlices
            % User played with rotations, slice order reversed (going down)
            % full K-space 3D encodings: descending
            seEncodings = linspace(1,-1,encodingPlan.numSE);
        else
            % slices are going up from 1 to the last
            % full K-space 3D encodings: ascending
            seEncodings = linspace(-1,1,encodingPlan.numSE);
        end
    else
        % single encoding, zero
        encodingPlan.numSE = 1;
        seEncodings = 0;
    end
else
    encodingPlan.numSE = 1;
    seEncodings = 0;
end

%% multiple echoes
encodingPlan.numCE          = encodingPlan.kSizeCE; % number of contrasts: sets of FE in each repetition
if contains(lower(acqData.pulseSeqFamilyName), 'gre-oop') || ...
        contains(lower(acqData.pulseSeqFamilyName), 'molli')
    encodingPlan.altMultiEcho   = 0; % Switch OFF Bipolar MultiEcho for GRE-OOP sequence
else
    encodingPlan.altMultiEcho   = 1; % multiple echoes in alternate direction
end
encodingPlan.encPerContrast = 1; % encodings (phase or slice) per contrast
encodingPlan.ceIncidence    = reshape(1:encodingPlan.numCE,[],1); % echo incidence

%% check parallel imaging or partial fourier acceleration
if ~strcmpi(acqData.parallelImaging, 'no')
    % parallel acceleration, for now only in the Phase encoding
    encodingPlan.rFactorPE = acqData.rFactor;
    encodingPlan.rFactorSE = 1;
    % no Fourier acceleration
    encodingPlan.fFactorFE = 1.0; % partial Fourier, <= 1.0
    encodingPlan.fFactorPE = 1.0; % partial Fourier, <= 1.0
    encodingPlan.fFactorSE = 1.0; % partial Fourier, <= 1.0
else
    % no parallel acceleration
    encodingPlan.rFactorPE = 1;
    encodingPlan.rFactorSE = 1;
    % check partial Fourier
    switch lower(acqData.partialFourier)
        case lower('readConjugate')
            % partial Fourier in the Frequency encoding (readout)
            encodingPlan.fFactorFE = acqData.fFactor; % partial Fourier, <= 1.0
            encodingPlan.fFactorPE = 1.0;
            encodingPlan.fFactorSE = 1.0; 
        case lower('phaseConjugate')
            % partial Fourier in the Phase encoding
            encodingPlan.fFactorFE = 1.0; 
            encodingPlan.fFactorPE = acqData.fFactor; % partial Fourier, <= 1.0
            encodingPlan.fFactorSE = 1.0; 
        case lower('sliceConjugate')
            % partial Fourier in the 3D encoding
            encodingPlan.fFactorFE = 1.0;
            encodingPlan.fFactorPE = 1.0;
            encodingPlan.fFactorSE = acqData.fFactor; % partial Fourier, <= 1.0
        otherwise
            encodingPlan.fFactorFE = 1.0; 
            encodingPlan.fFactorPE = 1.0;
            encodingPlan.fFactorSE = 1.0;
    end
end

%% index of starting encoding level depending on Fourier acc.
encodingPlan.conjFill = 0; % define if we are going to fill conj of K-space
% frequency encodings
numFE               = round(encodingPlan.fFactorFE*encodingPlan.numFE);
encodingPlan.startFE  = 1 + encodingPlan.numFE - numFE;
% phase encodings
numPE                   = round(encodingPlan.fFactorPE*encodingPlan.numPE);
encodingPlan.startPE    = 1 + encodingPlan.numPE - numPE;
% slice (3D) encodings
numSE                   = round(encodingPlan.fFactorSE*encodingPlan.numSE);
encodingPlan.startSE    = 1 + encodingPlan.numSE - numSE;

%% incidence arrays: give the position of acquired sample in the k-Space matrix
%   incidence of readouts in Frequency direction of K-space, useful for partial fourier
encodingPlan.feIncidence = ...
    reshape(encodingPlan.startFE:encodingPlan.numFE,[],1); 
% incidence of encoding in Phase direction of K-space: 
% include stride of parallel acc
encodingPlan.peIncidence = ...
    reshape(encodingPlan.startPE:encodingPlan.rFactorPE:encodingPlan.numPE,[],1); 
% incidence of encoding in Slice direction of K-space:
% include stride of parallel acc
encodingPlan.seIncidence = ...
    reshape(encodingPlan.startSE:encodingPlan.rFactorSE:encodingPlan.numSE,[],1);

%% select the encoding levels using the incidence matrices, from initial
% reduced encodings after acceleration
encodingPlan.peEncodings = reshape(peEncodings(encodingPlan.peIncidence),[],1);
encodingPlan.seEncodings = reshape(seEncodings(encodingPlan.seIncidence),[],1);

%% update number of encodings, including partial and parallel acc.
encodingPlan.numFE = length(encodingPlan.feIncidence);
encodingPlan.numPE = length(encodingPlan.peIncidence);
encodingPlan.numSE = length(encodingPlan.seIncidence);

%% acquisition info for reverse RO lines and PE incidence
encodingPlan.rxReverse   = zeros(encodingPlan.numPE*encodingPlan.numSE,1);
encodingPlan.rxReversePE = zeros(encodingPlan.numPE,1);
encodingPlan.rxReverseSE = zeros(encodingPlan.numSE,1);
if contains(lower(acqData.pulseSeqFamilyName), 'epi')
    encodingPlan.epiReverse             = 1; % reverse odd lines
    encodingPlan.rxReverse(1:2:end)     = 1; % TODO: for back compatibility, remove
    encodingPlan.rxReversePE(1:2:end)   = 1;
    encodingPlan.peIncidence            = flipud(encodingPlan.peIncidence);
else
    encodingPlan.epiReverse         = 0;
end

%% kSpace shift
if strcmpi(acqData.kspaceshift, 'yes')
    encodingPlan.kSpaceShift = 1;
else
    encodingPlan.kSpaceShift = 0;
end

%% fill conjugate K-space for SE sequences
if contains(lower(acqData.pulseSeqFamilyName), 'se')
    encodingPlan.conjFill = 1;
end

%% phase of the acquisition for each encoding 
encodingPlan.encPhase  = zeros(encodingPlan.numPE*encodingPlan.numSE,1); 
if contains(lower(acqData.pulseSeqFamilyName), 'bssfp') || ...
        strcmpi(acqData.pulseSeqFamilyName, 'molli')
    encodingPlan.encPhase(2:2:end) = pi; % change phase for BSSFP
end

if strcmpi(acqData.pulseSeqFamilyName, 'molli')
    encodingPlan.multipleSingleshot = 1;
else
    encodingPlan.multipleSingleshot = 0;
end

%% Normal encoding order: first phase then slice encodings
if contains(lower(acqData.pulseSeqFamilyName), 'rage') ... %MP-RAGE
        || contains(lower(acqData.pulseSeqFamilyName), 'conc') ... % concatenations
        || contains(lower(acqData.pulseSeqFamilyName), 'sh') % slice hopping
    encodingPlan.encPhaseDir = 'SP'; % change order, first slice, then phase
else
    encodingPlan.encPhaseDir = 'PS'; % first phase, then slice
end

%% reorder encodings for TSE sequences
encodingPlan.encPerShot = ones(numPE,1);
if contains(lower(acqData.pulseSeqFamilyName), 'tse') && (acqData.ETL > 1)
    numPE   = encodingPlan.numPE;
    numShot = ceil(numPE/acqData.ETL);
    encodingPlan.encPerShot = zeros(numShot,1);
    idxEnc  = zeros(numPE,1);
    counter = 0;
    for shot = 1:numShot
        % get the encoding indexes corresponding to this shot
        idxLocalEnc = shot:numShot:numPE;
        numLocalEnc = length(idxLocalEnc);
        % re-order shot encodings
        switch lower(acqData.encOrder)
            case 'centric'
                % order phase encoding and incidence in ABS ascending order
                [~,shotOrder] = sort( ...
                    abs(encodingPlan.peEncodings(idxLocalEnc)), 'ascend');
            case 'reversedcentric'
                % order phase encoding and incidence in ABS descencing order
                [~,shotOrder] = sort( ...
                    abs(encodingPlan.peEncodings(idxLocalEnc)), 'descend');
            case 'ascending'
                % order phase encoding and incidence in ascending order
                [~,shotOrder] = sort( ...
                    encodingPlan.peEncodings(idxLocalEnc), 'ascend');
            case 'descending'
                % order phase encoding and incidence in descencing order
                [~,shotOrder] = sort( ...
                    encodingPlan.peEncodings(idxLocalEnc), 'descend');
            case 'reversed'
                % reversed from default
                shotOrder = length(encodingPlan.peEncodings(idxLocalEnc)):-1:1; 
            otherwise
                % default: as is
                shotOrder = 1:length(encodingPlan.peEncodings(idxLocalEnc));
        end
        % store ordered encodings
        idxEnc(counter+1:counter+numLocalEnc) = idxLocalEnc(shotOrder);
        counter = counter + numLocalEnc;
        encodingPlan.encPerShot(shot) = numLocalEnc; % number of encodings in the shot
    end
    % apply the reorder to the encoding and incidence
    encodingPlan.peEncodings(:) = encodingPlan.peEncodings(idxEnc);
    encodingPlan.peIncidence(:) = encodingPlan.peIncidence(idxEnc);

else
    
    switch lower(acqData.encOrder)
        case 'centric'
            % order phase encoding and incidence in abs ascending order
            [~,idxEnc] = sort(abs(encodingPlan.peEncodings),'ascend');
        case 'reversedcentric'
            % order phase encoding and incidence in abs descencing order
            [~,idxEnc] = sort(abs(encodingPlan.peEncodings),'descend');
        case 'ascending'
            % order phase encoding and incidence in ascending order
            [~,idxEnc] = sort(encodingPlan.peEncodings,'ascend');   
        case 'descending'
            % order phase encoding and incidence in descending order
            [~,idxEnc] = sort(encodingPlan.peEncodings,'descend');
        case 'reversed'
            % reversed from default
            idxEnc = length(encodingPlan.peEncodings):-1:1;
        otherwise
            % default order: as is
            idxEnc = 1:length(encodingPlan.peEncodings);
    end
    % do not apply to EPI sequences
    if ~contains(lower(acqData.pulseSeqFamilyName), 'epi')
        encodingPlan.peEncodings(:) = encodingPlan.peEncodings(idxEnc);
        encodingPlan.peIncidence(:) = encodingPlan.peIncidence(idxEnc);
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
