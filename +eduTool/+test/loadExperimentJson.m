function [expControl,acquisition] = loadExperimentJson( inputSource )
%
% EDUTOOL.TEST.LOADEXPERIMENTJSON
%
%	Function that loads expControl and acquisition structs 
%   with data from an App-dependent source
%
% INPUT
%   inputSource   can be different types depending on the application
%   application   string with application type
%
% OUTPUT
%   expControl   expControl structure
%   acquisition  acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%


%% initialize the structure with defaults
[expControl]    = data.experiment.initializeExpControl();
[acquisition]   = data.experiment.initializeAcquisition();

%% load data from a valid json file inputSource
fid = fopen(inputSource,'r');
experimentData = jsondecode(fread(fid,inf,'*char').');
fclose(fid);

%% assign data by deep copy of experiment structures
% this will avoid crashes in case there are
% new fields that are not in the experiment data
expControl  = tools.misc.deepCopyStruct( ...
    expControl, experimentData.expControl);
acquisition = tools.misc.deepCopyStruct( ...
    acquisition, experimentData.acquisition);