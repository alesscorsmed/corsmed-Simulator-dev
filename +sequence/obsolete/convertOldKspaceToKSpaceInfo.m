function kSpaceInfo = convertOldKspaceToKSpaceInfo(info,...
    pulseSeqFamilyName,kSpaceInfo)

if nargin<2 || isempty(pulseSeqFamilyName)
    pulseSeqFamilyName = 'NA';
end

% incidence arrays: 
% give the position of acquired sample in the k-Space matrix
% incidence of readouts in FE direction of K-space, useful for partial fourier
kSpaceInfo.feIncidence  = transpose(1:kSpaceInfo.numFE);
% incidence of encoding in PE and SE direction of K-space + encoding order
if nargin<3
    %% encoding factors: size multipliers in the FE and PE directions
    kSpaceInfo.factorX = 1;
    kSpaceInfo.factorY = 1;

    %% K-space padding factors:
    %   if padded size is <= 1, no padding
    kSpaceInfo.xPadFactor = 1;
    kSpaceInfo.yPadFactor = 1;
    kSpaceInfo.zPadFactor = 1;

    % acceleration info
    kSpaceInfo.rFactorPE    = 1;
    kSpaceInfo.rFactorSE    = 1;
    kSpaceInfo.fFactorFE    = 1.0; % partial Fourier, <= 1.0
    kSpaceInfo.fFactorPE    = 1.0; % partial Fourier, <= 1.0
    kSpaceInfo.fFactorSE    = 1.0; % partial Fourier, <= 1.0
    
    kSpaceInfo.is3D  = 0;
    
    kSpaceInfo.peIncidence  = info.pulseSequence.kspace(:,6);
    kSpaceInfo.seIncidence  = info.pulseSequence.kspace(:,7);
    kSpaceInfo.encOrder     = 'PS'; % 'PS': phase-slice vs 'SP': slice-phase
    
    % acquisition info
    %   rxReverse: if 1, reverse RO (for EPI)
    kSpaceInfo.rxReverse    = zeros(kSpaceInfo.numPE*...
        kSpaceInfo.numSE*kSpaceInfo.numCN,1);
    %   shift for EPI
    kSpaceInfo.kSpaceShift  = 0;
    %   phase of the acquisition for each encoding
    kSpaceInfo.rxPhase      = info.pulseSequence.kspace(:,16);
    
    % start number for the encodings:
    %   starting number of the incidence, useful for Fourier acceleration
    kSpaceInfo.startFE      = 1;
    kSpaceInfo.startPE      = 1;
    kSpaceInfo.startSE      = 1;
    kSpaceInfo.startCE      = 1;
    
    %% final size of the image:
    kSpaceInfo.imSizeX = info.pulseSequence.Nx;
    kSpaceInfo.imSizeY = info.pulseSequence.Ny;
    kSpaceInfo.imSizeZ = max(info.pulseSequence.kspace(:,7));

    %% K-space sizes: usually original number of encodings x factors
    kSpaceInfo.xSize = info.pulseSequence.kspace(1,5);
    kSpaceInfo.ySize = max(info.pulseSequence.kspace(:,6));
    kSpaceInfo.zSize = max(info.pulseSequence.kspace(:,7));
    kSpaceInfo.cSize = max(info.pulseSequence.kspace(:,10));
    
    %% actual number of acquisition samples
    % this can be smaller than sizes above 
    % due to partial Fourier or parallel imaging
    kSpaceInfo.numFE        = info.pulseSequence.kspace(1,5);
    kSpaceInfo.numPE        = max(info.pulseSequence.kspace(:,6));
    kSpaceInfo.numSE        = max(info.pulseSequence.kspace(:,7));
    kSpaceInfo.numCN        = max(info.pulseSequence.kspace(:,10));   
    
else
    % avoid mismatch due to partial Fourier
    kSpaceInfo.peIncidence = transpose(1:size(kSpaceInfo.peIncidence));
    
    if strcmp(pulseSeqFamilyName,'SE') || strcmp(pulseSeqFamilyName,'IR-SE') || ...
        strcmp(pulseSeqFamilyName,'EPI') || strcmp(pulseSeqFamilyName,'SE-EPI')
            kSpaceInfo.peIncidence = flipud(kSpaceInfo.peIncidence);
    end
end
% incidence of encoding in Contrast direction of K-space
kSpaceInfo.coIncidence  = info.pulseSequence.kspace(:,10);