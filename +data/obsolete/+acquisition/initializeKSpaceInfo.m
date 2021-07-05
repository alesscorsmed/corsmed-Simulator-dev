function [kSpaceInfo] = initializeKSpaceInfo()
%
% DATA.ACQUISITION.INITIALIZEKSPACEINFO
%
%	Function that initializes the kSpace info structure.
%   Returns a kSpaceInfo with predefine values.
%   This structure has the info to assemble the K-space.
%
%   see reconstruction.signal.assembleKspace for how to use the data.
%
%
% INPUT
%   None
%
% OUTPUT
%   kSpaceInfo   kSpaceInfo structure
%
%========================  CORSMED AB © 2020 ==============================
%

%% final size of the image:
kSpaceInfo.imSizeX = 128;
kSpaceInfo.imSizeY = 128;
kSpaceInfo.imSizeZ = 1;

%% encoding factors: size multipliers in the FE and PE directions
kSpaceInfo.factorX = 1;
kSpaceInfo.factorY = 1;

%% K-space sizes: usually original number of encodings x factors
kSpaceInfo.xSize = 128;
kSpaceInfo.ySize = 128;
kSpaceInfo.zSize = 1;
kSpaceInfo.cSize = 1;
kSpaceInfo.is3D  = 0;

%% K-space padding factors:
%   if padded size is <= 1, no padding
kSpaceInfo.xPadFactor = 1;
kSpaceInfo.yPadFactor = 1;
kSpaceInfo.zPadFactor = 1;

%% actual number of acquisition samples
% this can be smaller than sizes above 
% due to partial Fourier or parallel imaging
kSpaceInfo.numFE = 128;
kSpaceInfo.numPE = 128;
kSpaceInfo.numSE = 1;
kSpaceInfo.numCN = 1; % number of contrasts: sets of FE in each repetition

% acceleration info
kSpaceInfo.rFactorPE = 1;
kSpaceInfo.rFactorSE = 1;
kSpaceInfo.fFactorFE = 1.0; % partial Fourier, <= 1.0
kSpaceInfo.fFactorPE = 1.0; % partial Fourier, <= 1.0
kSpaceInfo.fFactorSE = 1.0; % partial Fourier, <= 1.0

% start number for the encodings:
%   starting number of the incidence, useful for Fourier acceleration
kSpaceInfo.startFE = 1;
kSpaceInfo.startPE = 1;
kSpaceInfo.startSE = 1;

% incidence arrays: 
%   give the position of acquired sample in the k-Space matrix
% incidence of readouts in Frequency direction of K-space, useful for partial fourier
kSpaceInfo.feIncidence = zeros(kSpaceInfo.numFE,1);
 % incidence of encoding in Phase direction of K-space
kSpaceInfo.peIncidence = zeros(kSpaceInfo.numPE,1);
% incidence of encoding in Slice direction of K-space
kSpaceInfo.seIncidence = zeros(kSpaceInfo.numSE,1); 

% encodings
if kSpaceInfo.numPE == 1
    kSpaceInfo.peEncodings = 0;
else
    kSpaceInfo.peEncodings = linspace(-1,1,kSpaceInfo.numPE);
end
if kSpaceInfo.numSE == 1
    kSpaceInfo.seEncodings = 0;
else
    kSpaceInfo.seEncodings = linspace(-1,1,kSpaceInfo.numSE);
end

% acquisition info
%   rxReverse: if 1, reverse RO (for EPI)
kSpaceInfo.rxReverse = zeros(kSpaceInfo.numPE*kSpaceInfo.numSE,1);
%   shift for EPI
kSpaceInfo.kSpaceShift = 0;
%   phase of the acquisition for each encoding
kSpaceInfo.rxPhase   = zeros(kSpaceInfo.numPE*kSpaceInfo.numSE,1); 
%   encoding order
kSpaceInfo.encOrder = 'PS'; % 'PS': phase-slice vs 'SP': slice-phase
