function [expControl,acquisition] = fetchExperimentData(experiment,sessionData)
%
% EDUTOOL.SETUP.FETCHEXPERIMENTDATA
%
%	collects experiment data from DB.
%
% INPUT
%   experiment          struct with experiment info from frontend
%   sessionData         solution struct with initial data
%
% OUTPUT
%   expControl          initialized experiment control struct
%   acquisition         initialized acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.setup.fecthExperimentData';
if (nargin < 2)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% initialize expControl with sessionData
[expControl] = data.expControl.initialize(experiment,sessionData,'edutool');


%% Update ptx file based on parameters
% [expControl.simulation.kernelPtx] = eduTool.multiGPU.createPtxPath(...
%     expControl.model.pdInhomogeneity,expControl.motionSpecs.pattern,...
%      expControl.simulation.blocks,expControl.simulation.threads,...
%      expControl.simulation.simulationKernel);

% cuda architecture based on GPU
expControl.simulation.gpuName = gpuDevice(1).Name;
if contains(lower(expControl.simulation.gpuName),'v100')
    expControl.simulation.kernelPtx = sprintf('%s%s_sm70.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
end
if contains(lower(expControl.simulation.gpuName),'k80')
    expControl.simulation.kernelPtx = sprintf('%s%s_sm37.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
end

%% initialize acquisition data with sessionData
[acquisition] = data.acquisition.initialize(expControl,sessionData);

%% in case of motion or diffusion, disable time compression
if strcmpi(acquisition.data.pulseSeqFamilyName, 'pg-se-epi') ...
        || strcmpi(expControl.simulation.simulationEngine, 'diffusion') ...
        || ~strcmpi(expControl.motionSpecs.pattern, 'none')
    expControl.sequence.timeCompression = 0;
end

%% update related data
% update the number of simulations accordingly
if acquisition.data.is3D
    expControl.simulation.numberOfSim = 1;
else
    expControl.simulation.numberOfSim = acquisition.data.numSlices;
end
