function [ready,expControl,acquisition,inputSource] = fetchExperimentData(...
    inputSource,application,approach,mode)
%
% EDUTOOL.SETUP.FETCHEXPERIMENTDATA
%
%	collects experiment data from DB.
%
% INPUT
%   inputSource         source for the data: sessionData for eduTool
%   application         app type: 'eduTool' / 'apiJson' / 'expJson' / ...
%
% OUTPUT
%   expControl          initialized experiment control struct
%   acquisition         initialized acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.setup.fetchExperimentData';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

if nargin < 3
    approach    = '';
    mode        = '';
end    

%% load the experiment data using the corresponding App
[expControl,acquisition,inputSource] = data.experiment.loadExperiment(...
    inputSource,application,approach,mode);

%% check if experiment is ready
if isempty(expControl) || isempty(acquisition)
    % experiment not ready, return 0 flag
    ready = 0;
    return;
else
    ready = 1; % set flag as experiment ready
end

%% update acquisition with the computed Noise levels
[acquisition] = noise.calculateNoiseLevel( acquisition, expControl,...
    expControl.mrSystem.b0 );

%% in case we want to use old sequences, only analytical works
if strcmpi(acquisition.data.pulseSeqFamilyName, 'pg-se-epi')
    % Diffusion always uses new sequence
    expControl.useOldSequence = 0;
end
if ~isfield(expControl,'useOldSequence')
    expControl.useOldSequence = 0;
end
if expControl.useOldSequence
    expControl.simulation.simulationEngine = 'analytical';
    acquisition.data.encOrder              = 'std'; % no reorder availble
end

%% in case of motion or diffusion, disable time compression
if strcmpi(acquisition.data.pulseSeqFamilyName, 'pg-se-epi') ...
        || strcmpi(expControl.simulation.simulationEngine, 'diffusion') ...
        || ~strcmpi(expControl.motionSpecs.pattern, 'none')
    expControl.sequence.timeCompression = 0;
end

