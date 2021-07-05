function [simSignal] = engineV3Parallel( spinModel, pulseSequence, motionModel, ...
    expControl )
%
% EDUTOOL.RUN.ENGINEV3PARALLEL
%
%	Runs the simulator, uses parallel function evaluation
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.engineV3Parallel';
if (nargin < 4)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% if there is no DB connection
if ~isfield(expControl,'connLocalDB')
    expControl.connLocalDB = [];
end

%% SIMULATOR INSTANCE
expControl.estimatedRunTime = ...
    seconds(expControl.estimatedRunTime/expControl.simulation.numGPUs);
expControl.estimatedRunTime.Format = 'mm:ss';
eduToolMessage = sprintf('Pedal to the Metal... (Estimated Run Time %s sec)',...
    strrep(sprintf('%s',expControl.estimatedRunTime),':', ' min '));
eduTool.frontend.updateScannerStatus(expControl, eduToolMessage);

% intialize progress to 5%
expControl.progress = 5;

% remove conn and other data from expControl 
% to avoid issue with parallel pool feval
connLocalDB = expControl.connLocalDB;
expControl.connLocalDB = [];
gpuPool = expControl.simulation.gpuPool;
expControl.simulation.gpuPool = [];

%% extract simControl data from expControl
[simControl] = spinTwin.setup.initializeSimControl();
simControl.precision        = expControl.simulation.precision; % single / double
simControl.simulationEngine = expControl.simulation.simulationEngine; % Bloch / Phasor / PhasorPlus / Diffusion
simControl.odeMethod        = expControl.simulation.odeMethod; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
simControl.numGPUs          = expControl.simulation.numGPUs;
simControl.threads          = expControl.simulation.threads;
% correct analytical
if strcmpi(simControl.simulationEngine,'analytical')
    simControl.simulationEngine = 'Bloch';
    simControl.odeMethod        = 'analytical';
end

%% allocate output data (futures) for the solution of the simulation
fSim(1:spinModel.numJobs) = parallel.FevalFuture;
tSimulation(1:spinModel.numJobs) = tic();

%% evaluate the simulation in parallel
% loop on the number of Jobs
for jobNum = 1:spinModel.numJobs
    
    %% extract the model from the job queue
    model = spinModel.data{jobNum}.model;
    GPUindex = spinModel.data{jobNum}.GPUindex;
    currentGPU = gpuDevice;
    
    %% report start of simulation
    if expControl.debug.debugMode
        tSimulation(jobNum) = tic();
        fprintf(fid, ...
            '\n\n%s : starting simulation %d/%d, GPU index %d on GPU %d (%s)',...
            functionName, jobNum, spinModel.numJobs, GPUindex,...
            currentGPU.Index, currentGPU.Name);
        fprintf(fid, '\n');
    end
    
    %% run the simulation on the corresponding GPU
    fSim(jobNum) = parfeval(@spinTwin.runSimulation, 1, ...
        model, pulseSequence, motionModel, simControl);
    
end

%% prepare time array
rxTime = pulseSequence.time(pulseSequence.rxSignal(:) > 0);
% re-enable DB conn and other expControl data
expControl.connLocalDB = connLocalDB;
expControl.simulation.gpuPool = gpuPool;

%% loop on the number of Jobs
for ii = 1:spinModel.numJobs
    
    %% retrieve finished data when ready
    [jobNum, timeSolution] = fetchNext(fSim);

    %% store the solution
    simSignal.timeSolution{jobNum}.Sx   = timeSolution.Sx;
    simSignal.timeSolution{jobNum}.Sy   = timeSolution.Sy;
    simSignal.timeSolution{jobNum}.Sz   = timeSolution.Sz;
    simSignal.timeSolution{jobNum}.time = rxTime;
    
    %% store mapping info
    simSignal.timeSolution{jobNum}.sliceNum = spinModel.data{jobNum}.sliceNum;
    simSignal.timeSolution{jobNum}.numParts = spinModel.data{jobNum}.numParts;
    simSignal.timeSolution{jobNum}.partNum = spinModel.data{jobNum}.partNum;
    
    %% update progress bar (simulation assumed 90%)
    expControl.progress = expControl.progress +...
        90/spinModel.numJobs;
    eduTool.frontend.updateExperimentProgress(expControl);
    %% check status
    eduTool.frontend.checkExperimentStatus(...
        expControl,expControl.experimentID);
    
    %% report end of experiment
    if expControl.debug.debugMode
        fprintf(fid, ...
            '\n%s : done simulation %d/%d on GPU %d ',...
            functionName, jobNum, spinModel.numJobs,... 
            spinModel.data{jobNum}.GPUindex);
        fprintf(fid, '\n  Elapsed Time   %.3fs', toc(tSimulation(jobNum)));
        fprintf(fid, '\n');
    end
    
end

%% number of coils and reads
simSignal.numJobs   = spinModel.numJobs;
simSignal.numSlices = spinModel.numSlices;
simSignal.numCoils  = timeSolution.numCoils;
simSignal.numReads  = timeSolution.numReads;

%% update progress bar to 95%
expControl.progress = 95;
eduTool.frontend.updateExperimentProgress(expControl);

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Number of Jobs    %d', spinModel.numJobs);
    fprintf(fid, '\n  Number of Voxels  %d', spinModel.totalIso);
    fprintf(fid, '\n  Number of Steps   %d', pulseSequence.numSteps);
    fprintf(fid, '\n  Elapsed Time      %.3fs (%.1fps/voxel/step)',...
        tTotal, 1e12*tTotal/spinModel.totalIso/pulseSequence.numSteps);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end