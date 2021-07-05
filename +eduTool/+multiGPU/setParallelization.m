function [expControl,edtPool,gpuPool,numGPUs] = setParallelization( ...
    expControl,edtPool,gpuPool,numGPUs)
%
% EDUTOOL.MULTIGPU.SETPARALELLIZATION
%
%   updates the parallelization and kernels for the simulation.
%
%
% INPUT
%   expControl  struct with data
%   edtPool     current parPool struct
%   gpuPool     struct with gpuDevices
%   numGPUs     number of GPUs available
%
% OUTPUT
%   expControl  updated struct
%
% ========================= CORSMED AB @ 2020 ============================
%
 
%% make sure that we have correct parPool
if (nargin < 2) || isempty(edtPool)
    CPUperGPU = 1;
    myCluster = parcluster('local');
    numCPUs = myCluster.NumWorkers;
    numGPUs = gpuDeviceCount;
    if numGPUs > 0
        % initialize parPool if it does not exist
        numCPUs = min(numCPUs,CPUperGPU*numGPUs); % limit max of CPUs
        if isempty(gcp('nocreate'))
            edtPool = parpool(numCPUs, 'IdleTimeout', Inf);
        else
            edtPool = gcp();
        end
        % create a gpuPool with the GPU devices
        if edtPool.NumWorkers > 1
            spmd
                % each worker selects its own gpu
                gpuPool = gpuDevice(rem(labindex,numGPUs)+1);
            end
        else
            % single GPU, assign to current CPU
            gpuPool{1} = gpuDevice(labindex);
        end
    else
        MException('EduTool:BadInstance', 'Instance has not available GPUs');
    end
end
     
%% assign parPool to expControl
expControl.simulation.numCPUs = edtPool.NumWorkers;
expControl.simulation.numGPUs = numGPUs;
expControl.simulation.parPool = edtPool;
expControl.simulation.gpuPool = gpuPool;

%% get the GPU name to define the kernel to use
currentGPU = gpuPool{1};
expControl.simulation.gpuName = currentGPU.Name;

 %% cuda architecture based on GPU
% standard compiled with no flags
expControl.simulation.kernelPtx = sprintf('%s%s.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
% check GPU type and use optimized kernel
if contains(lower(expControl.simulation.gpuName),'v100')
    expControl.simulation.threads   = 256;
    expControl.simulation.kernelPtx = sprintf('%s%s_sm70.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
end
if contains(lower(expControl.simulation.gpuName),'k80')
    expControl.simulation.threads   = 128;
    expControl.simulation.kernelPtx = sprintf('%s%s_sm37.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
end
