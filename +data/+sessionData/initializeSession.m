function [sessionData] = initialize()
%
% DATA.SESSIONDATA.INITIALIZE
%
%	Function that initializes elements of the instance
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   Instance Attributes Structure
%
%========================  CORSMED AB © 2020 ==============================
%

% sessionData.anatModel.modelName               = '';
% sessionData.anatModel.defaultAnatModelPath    = '';
% sessionData.anatModel.defaultCoilsModelPath   = '';
% sessionData.anatModel.PDfluctuationFactor     = '';
% sessionData.anatModel.default_gridStep        = '';
% sessionData.anatModel.fatTissuesIDs           = '';
% sessionData.anatModel.PDinhomPathFile         = '';
% sessionData.anatModel.volumeGridSize          = '';

%% instance data
sessionData.versionNum     = 'v20200930a';
sessionData.instanceID     = '';
sessionData.instanceName   = '';
sessionData.courseID       = 0;
sessionData.userID         = 0;
sessionData.AWStagUserID   = 0;

%% drivers and dev info
sessionData.cudaVersion    = 8;
sessionData.parfeval       = 0;
sessionData.pythonVersion  = '3.7';
sessionData.developmentUse = 0;

%% folders
sessionData.folderSystem.testMode               = '';
sessionData.folderSystem.commonFolder           = '';
sessionData.folderSystem.kernelFolder           = '';
sessionData.folderSystem.resultsFolder          = '';
sessionData.folderSystem.gadgetronFolder        = '';
sessionData.folderSystem.gadgetronResultsFolder = '';
sessionData.folderSystem.gadgetronISMRMRDFolder = '';
sessionData.folderSystem.pulseSequenceFolder    = '';
sessionData.folderSystem.analyticReconFolder    = '';
sessionData.folderSystem.anatomicalModelFolder  = '';
sessionData.folderSystem.coilModelFolder        = '';
sessionData.folderSystem.errorFolder            = '';




