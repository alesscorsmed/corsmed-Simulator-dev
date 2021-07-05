function [encodingMap,encodingPlan] = oldKspaceToEncodingData(infoSequence,...
    pulseSeqFamilyName,encoding)

if nargin<2 || isempty(pulseSeqFamilyName)
    pulseSeqFamilyName = 'NA';
end

% incidence of encoding in PE and SE direction of K-space + encoding order
if nargin<3
    
    %% actual number of acquisition samples
    encodingPlan.numFE          = infoSequence.kspace(1,5);
    encodingPlan.numPE          = max(infoSequence.kspace(:,6));
    encodingPlan.numSE          = max(infoSequence.kspace(:,7));
    encodingPlan.numCE          = max(infoSequence.kspace(:,10));
    
    encodingPlan.kSizeFE        = infoSequence.kspace(1,5);
    encodingPlan.kSizePE        = max(infoSequence.kspace(:,6));
    encodingPlan.kSizeSE        = max(infoSequence.kspace(:,7));
    encodingPlan.kSizeCE        = max(infoSequence.kspace(:,10));
    
    encodingPlan.rFactorPE      = 1;
    encodingPlan.rFactorSE      = 1;
    
    % acquisition info
    %   rxReverse: if 1, reverse RO (for EPI)
    encodingPlan.rxReverse      = zeros(encodingPlan.numPE*...
        encodingPlan.numSE*encodingPlan.numCE,1);
    %   shift for EPI
    encodingPlan.kSpaceShift    = 0;
    
    % start number for the encodings:
    %   starting number of the incidence, useful for Fourier acceleration
    encodingPlan.startFE        = 1;
    encodingPlan.startPE        = 1;
    encodingPlan.startSE        = 1;
    
%     encodingPlan.feIncidence	= transpose(1:encodingPlan.numFE);
%     encodingPlan.peIncidence	= infoSequence.kspace(:,6);
%     encodingPlan.seIncidence	= infoSequence.kspace(:,7);
%     encodingPlan.ceIncidence    = infoSequence.kspace(:,10);
    
    encodingPlan.feIncidence    = ...
        reshape(encodingPlan.startFE:encodingPlan.numFE,[],1);
    encodingPlan.peIncidence    = ...
        reshape(encodingPlan.startPE:encodingPlan.rFactorPE:...
        encodingPlan.numPE,[],1);
    encodingPlan.seIncidence    = ...
        reshape(encodingPlan.startSE:encodingPlan.rFactorSE:...
        encodingPlan.numSE,[],1);    
    encodingPlan.ceIncidence    = ...
        reshape(1:encodingPlan.numCE,[],1); % echo incidence
    
    encodingPlan.encPerContrast = encodingPlan.numPE;
    
    encodingPlan.encPhaseDir	= 'PS'; % 'PS': phase-slice vs 'SP': slice-phase
    
    encodingPlan.encPhase       = infoSequence.kspace(:,16);
    
    encodingPlan.altMultiEcho   = 0;
    
    % Create a local expControl
    expControlLocal.debug.debugMode = 0;
    
    %% generate incidence matrices to map from time to k-space
    [encodingMap] = encoder.cartesian.encodingMapperKindex(...
        encodingPlan, expControlLocal);
    
%     %% final size of the image:
%     encodingMap.imSizeX     = infoSequence.Nx;
%     encodingMap.imSizeY     = infoSequence.Ny;
%     encodingMap.imSizeZ     = max(infoSequence.kspace(:,7));
%     
%     % acceleration info
%     encodingMap.rFactorPE 	= 1;
%     encodingMap.rFactorSE 	= 1;
%     encodingMap.fFactorFE	= 1.0; % partial Fourier, <= 1.0
	encodingPlan.factorX	= 1.0; % partial Fourier, <= 1.0
    encodingPlan.factorY	= 1.0; % partial Fourier, <= 1.0
    encodingPlan.fFactorZ 	= 1.0; % partial Fourier, <= 1.0
%     
    encodingPlan.is3D        = 0;
%     
%     %% encoding factors: size multipliers in the FE and PE directions
%     encodingMap.factorX     = 1;
%     encodingMap.factorY     = 1;
% 
    %% K-space padding factors:
    %   if padded size is <= 1, no padding
    encodingPlan.xPadFactor  = 1;
    encodingPlan.yPadFactor  = 1;
    encodingPlan.zPadFactor  = 1;

    
    
else
    % acquisition info
    %   rxReverse: if 1, reverse RO (for EPI)
    encoding.plan.rxReverse    = zeros(encoding.plan.numPE*...
        encoding.plan.numSE*encoding.plan.numCE,1);
    %   shift for EPI
    % incidence arrays: 
    % give the position of acquired sample in the k-Space matrix
    % incidence of readouts in FE direction of K-space, useful for partial fourier
    encoding.plan.feIncidence  = transpose(1:encoding.plan.numFE);
    % avoid mismatch due to partial Fourier
    encoding.plan.peIncidence = transpose(1:size(encoding.plan.peIncidence));
    
    if strcmp(pulseSeqFamilyName,'SE') || strcmp(pulseSeqFamilyName,'IR-SE') || ...
        strcmp(pulseSeqFamilyName,'EPI') || strcmp(pulseSeqFamilyName,'SE-EPI')
            encoding.plan.peIncidence = flipud(encoding.plan.peIncidence);
    end
    
    if encoding.plan.is3D
        encoding.plan.encPhaseDir = 'SP';
    end
    
    encodingPlan    = encoding.plan;
    encodingMap     = encoding.map;
    
end