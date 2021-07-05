function [expControl] = loadExpControl(expControl,experimentData,inputData)
%
% SERVICES.SETUP.LOADEXPCONTROL
%
%	Function that loads expControl struct from (REDIS) json
%
% INPUT
%   expControl   structure with acquisition parameters
%   experimentData   structure with json loaded parameters from Redis
%
% OUTPUT
%   expControl   structure with acquisition parameters
%
%========================  CORSMED AB © 2020 ==============================
%

if nargin<3
    inputData = [];
end

%% get info from the session
sessionData = experimentData.session;
expControl.versionNum       = 'S2V20210128';
expControl.userID           = sessionData.userId; 
expControl.instanceID       = sessionData.id; % assume this is instance?
expControl.courseID         = sessionData.labId; % now lab
expControl.skinID           = sessionData.skinId; % ## NOT in original expControl

expControl.anatomicalID     = inputData.anatomicalID;

if isfield(inputData,'approach') && strcmp(inputData.approach,'standalone')
    expControl.connLocalDB      = [];
end

%% Pass info for redis (if available)
if isfield(experimentData,'redis')
    expControl.redis = experimentData.redis;
end

%% debugging and development mode
% comes from CX, JV and GB instances
if ismember(expControl.userID,[790,933,1139])
    expControl.debug.devMode    = 1;
    expControl.debug.debugMode  = 1;
else
    expControl.debug.devMode    = 0;
    expControl.debug.debugMode  = 0;
end
expControl.debug.debugFile      = [];
% time stamps are created along with the allocation of expControl

%% folder system
% upgrade with more info
if isfield(inputData,'approach') && strcmp(inputData.approach,'jsonstandalone')
    if strcmp(inputData.mode,'dev')
        expControl.folderSystem.errorFolder         = sprintf('/home/ubuntu/edutoolTransferToS3/');
        expControl.folderSystem.statsFolder         = sprintf('/home/ubuntu/edutoolTransferToS3/');
        expControl.folderSystem.experimentFolder    = sprintf('/home/ubuntu/edutoolTransferToS3/');
        expControl.folderSystem.logFolder           = sprintf('/home/ubuntu/edutoolTransferToS3/');
    else
        expControl.folderSystem.errorFolder         = sprintf('/home/ubuntu/');
        expControl.folderSystem.statsFolder         = sprintf('/home/ubuntu/');
        expControl.folderSystem.experimentFolder    = sprintf('/home/ubuntu/');
        expControl.folderSystem.logFolder           = sprintf('/home/ubuntu/');
    end
    expControl.folderSystem.highPriorityFolder      = sprintf('/home/ubuntu/reconimages/');
else
    mainLocalFolder     = '/efs-mount-point-MATLAB/cxanthis/TEMP';
    mainRemoteFolder    = '/efs-mount-point-MATLAB/cxanthis/TEMP';
    expControl.folderSystem.baseFolder = mainLocalFolder;
    % update kernel path
    expControl.simulation.kernelFolder = ...
        sprintf('%s/PROJECTS/%s/kernels/',...
        expControl.folderSystem.baseFolder, expControl.application);
    % default kernel
    expControl.simulation.kernelPtx = sprintf('%s%s.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
    % upgrade with more info
    expControl.folderSystem.errorFolder = ...
        sprintf('%s/ERRORS/%s/', ...
        expControl.folderSystem.baseFolder, expControl.application);
    expControl.folderSystem.experimentFolder = ...
        sprintf('%s/RESULTS/%s/user_%d/', ...
        expControl.folderSystem.baseFolder, ...
        expControl.application, expControl.userID);
    expControl.folderSystem.logFolder = sprintf('%s/LOGS/%s/', ...
        expControl.folderSystem.baseFolder, expControl.application);
    if ~isfolder(expControl.folderSystem.experimentFolder)
        mkdir(expControl.folderSystem.experimentFolder);
    end
end

%% get info from experiment
experiment = experimentData.experiment;
expControl.experimentID     = experiment.id;
expControl.remotedbID       = [];  % ## NOT available in experiment
expControl.status           = experiment.status;
expControl.experInfo        = experiment.info;
expControl.reconstructor    = experiment.reconstructor;
expControl.reconInfo        = experiment.reconInfo;
expControl.pulseqID         = experiment.pulseSequenceID;
expControl.pathPulseSeq     = ''; % NOTE: this was generated, but it is useless in modular 
expControl.pulseqMatFile    = experiment.pulseqMatFileName;

expControl.other.pulseqMatFileID  = experiment.pulseqMatFileId; % ## NOT in original expControl
expControl.other.resultsPath      = experiment.resultsPath; % ## NOT in original expControl
expControl.other.reconPath        = experiment.reconPath; % ## NOT in original expControl
expControl.other.sar              = experiment.sar; % ## NOT in original expControl
expControl.other.sarBackend       = experiment.sarBackend; % ## NOT in original expControl
expControl.other.configurations   = experiment.configurations; % ## NOT in original expControl
expControl.other.execTimeSec      = experiment.execTimeSec; % ## NOT in original expControl
expControl.other.execInInstance   = experiment.execInInstance; % ## NOT in original expControl

%% parameters of the configuration
dataConfiguration = experimentData.configuration;
for i=1:length(dataConfiguration)
    structExper.(dataConfiguration(i).name) = dataConfiguration(i).selectedValue;
end
% assign
expControl.outerFOVratio = str2double(structExper.outer_fov_ratio);

% general
expControl.edutoolVersion           = structExper.edutool_version;
expControl.sequenceVersion          = structExper.pulse_seq_generator;
expControl.textOnViewers            = str2double(structExper.text_on_viewers);
expControl.advancedNotifications    = str2double(structExper.advanced_notifications);

% simulation related
expControl.simulation.analyticalSim     = str2double(structExper.analytical_simulation); % full analytical
expControl.simulation.simulationEngine  = structExper.simulationMode; % numerical / analytical
expControl.simulation.simulationKernel  = structExper.simulationKernel; % latest / stable
expControl.simulation.threads           = str2double(structExper.kernel_threads);
expControl.simulation.blocks            = str2double(structExper.kernel_blocks);
expControl.simulation.activateCS        = str2double(structExper.activate_cs);
expControl.simulation.applyNoise        = str2double(structExper.noise);
expControl.simulation.testSpurious      = str2double(structExper.testSpuriousEchoes);
expControl.simulation.activateSusc      = str2double(structExper.susceptibility);

% model related
expControl.model.selectedRegion     = str2double(structExper.selected_region);
expControl.model.zIsocenter         = str2double(structExper.isocenter)/100;
if (expControl.courseID == 6)
    % MIDA, fix isocenter, need to correct in Front End
    expControl.model.zIsocenter = 0.12;
end
expControl.model.zIsocenterLimit    = str2double(structExper.isocenter_limit);
expControl.model.coilType           = structExper.coil;
expControl.model.coilMode           = structExper.coil_basic_or_advanced;  % It refers to the option coil_basic_or_advanced
expControl.model.isotropicGrid      = str2double(structExper.isotropic_grid);
expControl.model.gridSizeOption     = structExper.gridSizeOption;
switch lower(expControl.model.gridSizeOption)
    case '1mm'
        expControl.model.gridStep   = [1.0, 1.0, 1.0]*1e-3;
    case '0p5mm'
        expControl.model.gridStep   = [0.5, 0.5, 0.5]*1e-3;
    otherwise
        expControl.model.gridSizeOption = 'optimal';
        expControl.model.gridStep   = [0.5, 0.5, 1.0]*1e-3;
        expControl.simulation.simulationEngine = 'phasor';
end
% for native model resolution, use phasor
if ~expControl.model.isotropicGrid
    expControl.simulation.simulationEngine = 'phasor';
end
expControl.model.useSliceThickness  = str2double(structExper.performance_gridZ_sliceThickness);
expControl.model.pdInhomogeneity    = str2double(structExper.PDinhomogeneity);


% sequence related
expControl.sequence.timeCompression = str2double(structExper.fast_algorithm); % former fast algorithm
expControl.sequence.dwellTime       = structExper.dwell_time;
expControl.sequence.deactivateGx    = str2double(structExper.deactivateGx);
expControl.sequence.deactivateGy    = str2double(structExper.deactivateGy);
expControl.sequence.deactivateGz    = str2double(structExper.deactivateGz);
% Use the old (v1) or new (v2) pulse sequence generator
pulseSeqGenerator = structExper.pulse_seq_generator;
if strcmpi(pulseSeqGenerator,'v1')
    expControl.useOldSequence = 1;
else
    expControl.useOldSequence = 0;
end

% motion related
expControl.motionSpecs.pattern      = structExper.motion;
expControl.motionSpecs.rotAngle     = str2double(structExper.motion_rot_angle);
expControl.motionSpecs.rotFreq      = str2double(structExper.motion_rot_freq);
expControl.motionSpecs.transMag     = str2double(structExper.motion_trans_magn);
expControl.motionSpecs.transFreq    = str2double(structExper.motion_trans_freq);
expControl.motionSpecs.transAxis    = str2double(structExper.motion_trans_axis);

%% TEST
% If anatomicalID = 9, get from db the default cardiac phase
if (expControl.anatomicalID == 9)
    expControl.model.cardiacPhase  = str2double(structExper.defaultCardiacPhase);
end

%%


% standard and other
expControl.advancedNotifications    = str2double(structExper.advanced_notifications);

% non-existent before
expControl.other.anatomical             = str2double(structExper.anatomical);  
expControl.other.monitorPerformance     = str2double(structExper.monitor_performance);
expControl.other.centralDB              = str2double(structExper.central_db);
expControl.other.openStackModel         = str2double(structExper.open_stack_mode);


