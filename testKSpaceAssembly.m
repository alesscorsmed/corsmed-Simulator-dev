
%% modify this info to see the effect
expControl.debug.debugMode = 1;

sizeX = 128;
sizeY = 60;
sizeZ = 3;
sizeC = 10;

rFactor         = 1;
fourierFactorX  = 1;
fourierFactorY  = 1;
fourierFactorZ  = 1.0;

typeBSSFP           = 0; % alternate phase
typeEPI             = 0; % reverses readout direction each encoding
altMultiEcho        = 0; % reverses readout direction each echo
numEncPerContrast   = 180; % number of phase/slice encodings per contrast
encPhaseDir         = 'PS';

%% create the relevant encoding plan info
encodingPlan.kSizeFE = sizeX;
encodingPlan.kSizePE = sizeY;
encodingPlan.kSizeSE = sizeZ;
encodingPlan.kSizeCE = sizeC;

%% incidences
encodingPlan.numFE    = round(fourierFactorX*sizeX);
encodingPlan.startFE  = 1 + sizeX - encodingPlan.numFE;
encodingPlan.numPE    = round(fourierFactorY*sizeY);
encodingPlan.startPE  = 1 + sizeY - encodingPlan.numPE;
encodingPlan.numSE    = round(fourierFactorZ*sizeZ);
encodingPlan.startSE  = 1 + sizeZ - encodingPlan.numSE;
encodingPlan.numCE    = sizeC;


%% incidence arrays: give the position of acquired sample in the k-Space matrix
encodingPlan.feIncidence = ...
    reshape(encodingPlan.startFE:sizeX,[],1); 
encodingPlan.peIncidence = ...
    reshape(encodingPlan.startPE:rFactor:sizeY,[],1); 
encodingPlan.seIncidence = ...
    reshape(encodingPlan.startSE:sizeZ,[],1);
encodingPlan.ceIncidence = ...
    reshape(1:sizeC,[],1); 

% reverse the readout direction
rxReverse = zeros(encodingPlan.numPE*encodingPlan.numSE,1);
rxPhase   = zeros(encodingPlan.numPE*encodingPlan.numSE,1);

%% BSSFP case
if typeBSSFP
    % tx/rx phase
    rxPhase(2:2:end) = pi;
end

%% EPI case
if typeEPI
    kSpaceShift = 1;
    % reverse readout direction
    rxReverse(2:2:end) = 1;
else
    kSpaceShift = 0;
end

%% multi contrast
% extend rx arrays to number of contrast
rxReverse   = repmat(rxReverse, 1, encodingPlan.numCE );
rxPhase     = repmat(rxPhase,   1, encodingPlan.numCE );

% modify to multi echo
if altMultiEcho
    rxReverse(:,2:2:end) = ~rxReverse(:,2:2:end);
end

% reorganize depending on the number of encodings per contrast
rxReverse   = reshape(rxReverse, numEncPerContrast, [], encodingPlan.numCE);
rxReverse   = permute(rxReverse, [1,3,2]);
rxPhase     = reshape(rxPhase, numEncPerContrast, [], encodingPlan.numCE);
rxPhase     = permute(rxPhase, [1,3,2]);

%% assign data to encoding plan
encodingPlan.kSpaceShift	= kSpaceShift;
encodingPlan.encPhaseDir	= encPhaseDir;
encodingPlan.encPerContrast = numEncPerContrast;
encodingPlan.rxReverse      = reshape(rxReverse,[],1);
encodingPlan.encPhase   	= reshape(rxPhase,[],1);

%% Return the mapper
[encodingMap] = encoder.cartesian.encodingMapperKindex( encodingPlan, expControl );

kSpaceInc = reshape(encodingMap.kSpaceIdx, ...
    encodingMap.kSizeFE, encodingMap.kSizePE, encodingMap.kSizeSE, encodingMap.kSizeCE);

