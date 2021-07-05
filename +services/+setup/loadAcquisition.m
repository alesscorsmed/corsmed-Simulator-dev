function [acquisition,experimentData] = ...
    loadAcquisition(acquisition,experimentData)
%
% SERVICES.SETUP.LOADACQUISITION
%
%	Function that loads acquisition struct from (REDIS) json
%
% INPUT
%   acquisition   structure with acquisition parameters
%   experimentData   structure with json loaded parameters from Redis
%
% OUTPUT
%   acquisition   structure with acquisition parameters
%
%========================  CORSMED AB © 2020 ==============================
%


%% pulse sequence data
dataPulseSequence = experimentData.pulseSequence;
pulseSeqSchem       = dataPulseSequence.pulseSequenceSchemeId;
pulseSeqNum         = dataPulseSequence.seqnum;
pointsAll           = dataPulseSequence.simPoints;
pointsAllFrontend   = dataPulseSequence.points;
%% NOTE: need test for multiple slices, how the data will come
if ~isempty(pointsAll)
    pointsAll           = strsplit(pointsAll,'|');
end
if ~isempty(pointsAllFrontend)
    pointsAllFrontend   = strsplit(pointsAllFrontend,'|');
end
% extra data:
% dataPulseSequence.seqnum: []
% dataPulseSequence.pulseSequenceSchemeId: 1
% dataPulseSequence.id: 3
% dataPulseSequence.labId: 1
% dataPulseSequence.userID: 1
% dataPulseSequence.clonedPulseqId: []
% dataPulseSequence.info: []
% dataPulseSequence.name: 'test'
% dataPulseSequence.points: []
% dataPulseSequence.irlInfo: []
% dataPulseSequence.snr: 0
% dataPulseSequence.noiseStd: 0
% dataPulseSequence.batchId: []
% dataPulseSequence.simPoints: []
% dataPulseSequence.cameFrom: []
% dataPulseSequence.createdInVersion: []
% dataPulseSequence.folder: []
% dataPulseSequence.liveIrl: []
% dataPulseSequence.liveSar: []

%% pulse sequence scheme
dataPulseScheme = experimentData.pulseSequenceScheme;
pulseSeqType        = dataPulseScheme.typeId;
pulseSeqFamilyName  = dataPulseScheme.filename;
pulseSeqName        = dataPulseScheme.name;
pulseSeqAnalytical  = dataPulseScheme.analytical;

% Tread the GRE-concatenations as a 3D sequence so as the simulator to
% treat the anatomical model of interest as a slab instead of a slice (for
% interpolation purposes)
if contains(lower(pulseSeqFamilyName), 'conc')
    is3D = 1;
    experimentData.pulseSequenceScheme.is3D = is3D;
else
    is3D = experimentData.pulseSequenceScheme.is3D;
end

% other data:
% dataPulseScheme.id: 1
% dataPulseScheme.typeId: 2
% dataPulseScheme.active: 1
% dataPulseScheme.filename: 'GRE'
% dataPulseScheme.name: 'Gradient Echo (GRE)'
% dataPulseScheme.description: '<h3>Gradient Echo (GRE)</h3><p>Gradient echo (GRE) is a class of pulse sequences that is primarily used for fast scanning. GRE is widely used in applications that require acquisition speed. This is accomplished due to the low flip angle of the excitation pulse, which is typically less than 90. Therefore, GRE pulse sequences can use short TR (2-50 ms).</p>'
% dataPulseScheme.analytical: 1
% dataPulseScheme.is3D: 0

%% parameters of the sequence
dataSequenceValues = experimentData.pulseSequenceValues;
for i=1:length(dataSequenceValues)
    structExper.(dataSequenceValues(i).uniqueName) = dataSequenceValues(i).selectedValue;
end

%% noise
dataNoiseLevels = experimentData.noiseLevels;
refVoxelSizeX   = str2double(dataNoiseLevels.voxelsizexref);
refVoxelSizeY   = str2double(dataNoiseLevels.voxelsizeyref);
refVoxelSizeZ   = str2double(dataNoiseLevels.voxelsizezref);
refBW           = str2double(dataNoiseLevels.bwref);
refB0           = str2double(dataNoiseLevels.b0ref);
refNEX          = dataNoiseLevels.nexref;
refSizeX        = dataNoiseLevels.kspacexref;
refSizeY        = dataNoiseLevels.kspaceyref;
refSizeZ        = dataNoiseLevels.kspacezref; % Number of encoding slices
noiseLevel      = str2double(dataNoiseLevels.noiseLevel); % actual noise level
if isnan(noiseLevel)
    noiseLevel = 0.0;
end
refVoxelVolume  = refVoxelSizeX*refVoxelSizeY*refVoxelSizeZ;
refEncSize      = refSizeX*refSizeY*refSizeZ;
% other data:
% dataNoiseLevels.id: 1
% dataNoiseLevels.pulseSequenceSchemeID: 1
% dataNoiseLevels.voxelsizexref: '0.001'
% dataNoiseLevels.voxelsizeyref: '0.001'
% dataNoiseLevels.voxelsizezref: '0.001'
% dataNoiseLevels.bwref: '200000'
% dataNoiseLevels.b0ref: '1.5'
% dataNoiseLevels.nexref: 1
% dataNoiseLevels.kspacexref: 100
% dataNoiseLevels.kspaceyref: 100
% dataNoiseLevels.kspacezref: 0
% dataNoiseLevels.noiseLevel: ''
    
%% process the data and incorporate into the acquisition struct

% sequence data
acquisition.data.pulseSeqSchem       = pulseSeqSchem;
acquisition.data.pulseSeqNum         = pulseSeqNum;
acquisition.data.pulseSeqName        = pulseSeqName;
acquisition.data.pulseSeqFamilyName  = pulseSeqFamilyName;
acquisition.data.pulseSeqType        = pulseSeqType;
acquisition.data.pulseSeqAnalytical  = pulseSeqAnalytical;
acquisition.data.is3D                = is3D;

if strcmp(pulseSeqName,'DIXON')
    acquisition.data.pulseSeqActualName = pulseSeqName;
else
    acquisition.data.pulseSeqActualName = pulseSeqFamilyName;
end

% plane
acquisition.data.pointsAll           = pointsAll;
acquisition.data.pointsAllFrontend   = pointsAllFrontend;

% basic info, should always be available
acquisition.data.numSlices      = str2num(structExper.slices);
acquisition.data.sliceThickness = str2num(structExper.sliceThickness);
if is3D && ~contains(lower(pulseSeqFamilyName), 'conc')
    acquisition.data.sliceGap = 0.0;
    acquisition.data.fovSE = acquisition.data.numSlices *...
        ( acquisition.data.sliceThickness + acquisition.data.sliceGap );
    acquisition.data.numSE = acquisition.data.numSlices;
    acquisition.data.matrixZ = acquisition.data.numSlices;
else
    acquisition.data.sliceGap = str2num(structExper.slice_gap);
    acquisition.data.fovSE = acquisition.data.sliceThickness;
    acquisition.data.numSE = 0;
    acquisition.data.matrixZ = 1;
end

acquisition.data.matrixX = str2num(structExper.matrix_x);
acquisition.data.matrixY = str2num(structExper.matrix_y);

acquisition.data.fovFE = str2num(structExper.fov_x);
acquisition.data.fovPE = str2num(structExper.fov_y);

acquisition.data.numFE = acquisition.data.matrixX;
acquisition.data.numPE = acquisition.data.matrixY;

if isfield(structExper, 'te')
    acquisition.data.TE = str2num(structExper.te);
end

if isfield(structExper, 'tr')
    acquisition.data.TR = str2num(structExper.tr);
end
if isfield(structExper, 'mp_rage_tr')
    acquisition.data.shotTR = str2num(structExper.mp_rage_tr);
else
    acquisition.data.shotTR = acquisition.data.TR;
end

acquisition.data.rxBW    = str2num(structExper.BW);
if isfield(structExper, 'BWpx')
    acquisition.data.pixelBW = str2num(structExper.BWpx);
end

acquisition.data.NEX = str2num(structExper.NEX);

acquisition.data.reconstructor = structExper.reconstructor;

% main RF
if isfield(structExper, 'RFflipAngle')
    acquisition.mainRF.flipAngle       = str2num(structExper.RFflipAngle);
else
    acquisition.mainRF.flipAngle = 90;
end
acquisition.mainRF.duration        = str2num(structExper.RFduration);
acquisition.mainRF.sliceThickness  = str2num(structExper.sliceThickness);
% refocusing RF consistent if needed
acquisition.refRF.duration        = str2num(structExper.RFduration);
acquisition.refRF.sliceThickness  = str2num(structExper.sliceThickness);

% optional entries
if isfield(structExper, 'OOP_deltaTime')
    if (isfield(structExper,'b0map') && strcmp(structExper.b0map,'1'))
        acquisition.data.numEchoes = 3;
        acquisition.data.generateB0map = 1;
    else
        acquisition.data.numEchoes = 2;
        acquisition.data.generateB0map = 0;
    end
    acquisition.data.deltaTimeOOP = str2num(structExper.OOP_deltaTime);
end
if strcmpi(pulseSeqFamilyName, 'molli')
    acquisition.data.MOLLI.scheme   = structExper.molli_scheme;
    acquisition.data.MOLLI.TIs      = structExper.tis;
    
    schemeCellArray = strsplit(acquisition.data.MOLLI.scheme,'(');
    LLexper = 0;
    schemeMOLLI = zeros(1,size(schemeCellArray,2));
    for iLL = 1:size(schemeCellArray,2)
        LLexper = LLexper + 1;
        schemeMOLLI(1,LLexper)      = str2double(schemeCellArray{1,iLL}(end));
    end
    acquisition.data.numEchoes      = sum(schemeMOLLI);
end

if isfield(structExper, 'kspaceshift')
    acquisition.data.kspaceshift  = structExper.kspaceshift;
end
if isfield(structExper, 'foldoversuppr')
    acquisition.data.foldoverSuppr  = structExper.foldoversuppr;
    if strcmpi(acquisition.data.foldoverSuppr, 'yes')
        acquisition.data.samplingFactorPE   = 2;
    else
        acquisition.data.samplingFactorPE   = 1;
    end
end
if isfield(structExper, 'foldoverdir')
    acquisition.data.foldoverDir    = structExper.foldoverdir;
end
if isfield(structExper, 'etl')
    acquisition.data.ETL = str2num(structExper.etl);
end
if isfield(structExper, 'concatenations')
    acquisition.data.concatenations = str2num(structExper.concatenations);
end
if isfield(structExper, 'parallelImaging')
    acquisition.data.parallelImaging = structExper.parallelImaging;
    if isfield(structExper, 'rfactor')
        acquisition.data.rFactor = str2num(structExper.rfactor);
    end
    if strcmpi(structExper.parallelImaging, 'sense')
        acquisition.data.reconstructor   = 'sense';
    end
end
if isfield(structExper, 'partialFourier')
    acquisition.data.partialFourier = structExper.partialFourier;
end
if isfield(structExper, 'partialFourierFactor')
    acquisition.data.fFactor = str2num(structExper.partialFourierFactor);
end
if isfield(structExper, 'encodingOrder')
    acquisition.data.encOrder = structExper.encodingOrder;
end

if isfield(structExper, 'refoc_pulse_angle')
    acquisition.refRF.flipAngle = str2num(structExper.refoc_pulse_angle);
end
if isfield(structExper, 'ir_time')
    acquisition.prepIR.Apply     = 1;
    acquisition.prepIR.TI        = str2num(structExper.ir_time);
end
if isfield(structExper, 'preparation_rf_angle')
    acquisition.prepIR.flipAngle = str2num(structExper.preparation_rf_angle);
else
    acquisition.prepIR.flipAngle = 180; % in degrees
end
if isfield(structExper, 'ir_duration')
    acquisition.prepIR.duration  = str2num(structExper.ir_duration);
end
if isfield(structExper, 'ir_type')
    acquisition.prepIR.type      = structExper.ir_type;
end
if isfield(structExper, 'fatsat')
    acquisition.data.fatsat  = structExper.fatsat;
end
% PG
if isfield(structExper, 'PG_AG')
    acquisition.encPG.AG  = str2num(structExper.PG_AG)*1e-3;
end
if isfield(structExper, 'PG_TG')
    acquisition.encPG.TG  = str2num(structExper.PG_TG);
end
if isfield(structExper, 'PG_TG')
    acquisition.encPG.Dir  = structExper.PG_DG;
end
acquisition.encPG.TAU  = acquisition.data.TE;

if isfield(structExper, 'vps')
    acquisition.data.vps = str2num(structExper.vps);
end

% noise data
acquisition.noise.refVoxelVolume= refVoxelVolume;
acquisition.noise.refEncSize    = refEncSize;
acquisition.noise.refNEX        = refNEX;
acquisition.noise.refBW         = refBW;
acquisition.noise.refB0         = refB0;
acquisition.noise.noiseLevel    = noiseLevel;

