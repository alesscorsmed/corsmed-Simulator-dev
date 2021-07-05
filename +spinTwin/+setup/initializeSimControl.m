function [simControl] = initializeSimControl()
%
% Initialize the simControl struct
%
%========================  CORSMED AB © 2020 ==============================

%% load json file with version info
kernelVersionFile = '+spinTwin/kernelVersion.json';
fid = fopen(kernelVersionFile,'r');
kernelVersion = jsondecode(fread(fid,inf,'*char').');
fclose(fid);
% add the version and kernels
simControl.version          = kernelVersion.version;   % current version
simControl.refKernel        = kernelVersion.refKernel; % stable verified kernel
simControl.simKernel        = kernelVersion.simKernel; % kernel for forward simulation
simControl.sbrKernel        = kernelVersion.sbrKernel; % kernel for sbr simulation
simControl.devKernel        = kernelVersion.devKernel; % kernel for dev projects (wildcard useful for develop/test/compare)
simControl.kernelPath       = kernelVersion.kernelPath;

%% basic info
simControl.precision        = 'single'; % single / double
simControl.simulationEngine = 'Bloch'; % Bloch / Phasor / Diffusion
simControl.odeMethod        = 'analytical'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
simControl.numGPUs          = 1;
simControl.threads          = 256;

%% cuda architecture and ptx path based on GPU
currentGPU = gpuDevice();
% check GPU type and use optimized kernel
if contains(lower(currentGPU.Name),'v100')
    simControl.ptxSIM = sprintf('%s%s_sm70.ptx', ...
        simControl.kernelPath, simControl.simKernel);
    simControl.ptxSBR = sprintf('%s%s_sm70.ptx', ...
        simControl.kernelPath, simControl.sbrKernel);
    simControl.ptxDEV = sprintf('%s%s_sm70.ptx', ...
        simControl.kernelPath, simControl.devKernel);
end
if contains(lower(currentGPU.Name),'k80')
    simControl.ptxSIM = sprintf('%s%s_sm37.ptx', ...
        simControl.kernelPath, simControl.simKernel);
    simControl.ptxSBR = sprintf('%s%s_sm37.ptx', ...
        simControl.kernelPath, simControl.sbrKernel);
    simControl.ptxDEV = sprintf('%s%s_sm37.ptx', ...
        simControl.kernelPath, simControl.devKernel);
end
