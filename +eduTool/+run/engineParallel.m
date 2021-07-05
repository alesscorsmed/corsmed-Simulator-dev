function [simSignal] = engineParallel( spinModel, pulseSequence, expControl )
%
% EDUTOOL.RUN.ENGINEPARALLEL
%
%	Runs the simulator, uses parallel function evaluation
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.engineParallel';
if (nargin < 3)
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

%% prepare time array
rxTime = pulseSequence.time(pulseSequence.rxSignal(:) > 0);
    
%% double loop approach to limit parallel execution and memory
numGPUs     = expControl.simulation.numGPUs;
numCPUs     = expControl.simulation.numCPUs;
numPROC     = max(numGPUs,numCPUs);
numSEQIT    = ceil(spinModel.numJobs/numPROC); % number of sequential iterations

%% external sequential loop
for kk = 1:numSEQIT
    
    % remove conn and other data from expControl
    % to avoid issue with parallel pool feval
    connLocalDB = expControl.connLocalDB;
    expControl.connLocalDB = [];
    gpuPool = expControl.simulation.gpuPool;
    expControl.simulation.gpuPool = [];

    %% allocate output data (futures) for the solution of the simulation
    localNumProc = min(numPROC, spinModel.numJobs - numPROC*(kk-1)); % remaining jobs
    fSim(1:localNumProc) = parallel.FevalFuture;
    tSimulation(1:localNumProc) = tic();
    if expControl.debug.debugMode
        fprintf(fid, '\n\n');
    end
    
    %% evaluate the simulation in parallel
    % loop on the number of parallel processes
    for pp = 1:localNumProc
        
        %% extract the model from the job queue
        jobNum = pp + numPROC*(kk-1);
        model = spinModel.data{jobNum}.model;
        GPUindex = spinModel.data{jobNum}.GPUindex;
        currentGPU = gpuDevice;
        
        %% report start of simulation
        if expControl.debug.debugMode
            tSimulation(pp) = tic();
            fprintf(fid, ...
                '%s : (%s) starting simulation %d/%d, GPU index %d on GPU %d (%s)',...
                functionName, expControl.simulation.simulationEngine,...
                jobNum, spinModel.numJobs, GPUindex,...
                currentGPU.Index, currentGPU.Name);
            fprintf(fid, '\n');
        end
        
        %% run the simulation on the corresponding GPU
        fSim(pp) = parfeval(@simulator.runSimulation, 1, ...
            model, pulseSequence, expControl, [], GPUindex);
        %% check status
        eduTool.frontend.checkExperimentStatus(...
            expControl,expControl.experimentID);
    end
    
    % re-enable DB conn and other expControl data
    expControl.connLocalDB = connLocalDB;
    expControl.simulation.gpuPool = gpuPool;
    if expControl.debug.debugMode
        fprintf(fid, '\n');
    end
    
    %% loop on the number of parallel Jobs
    for pp = 1:localNumProc
                
        %% retrieve finished data when ready
        [jobId, timeSolution] = fetchNext(fSim);
        
        %% store the solution
        jobNum = jobId + numPROC*(kk-1);
        simSignal.timeSolution{jobNum}.Sx   = timeSolution.Sx;
        simSignal.timeSolution{jobNum}.Sy   = timeSolution.Sy;
        simSignal.timeSolution{jobNum}.Sz   = timeSolution.Sz;
        simSignal.timeSolution{jobNum}.time = rxTime;
        
        %% store mapping info
        simSignal.timeSolution{jobNum}.phaseNum = spinModel.data{jobNum}.phaseNum;
        simSignal.timeSolution{jobNum}.contrNum = spinModel.data{jobNum}.contrNum;
        simSignal.timeSolution{jobNum}.frameNum = spinModel.data{jobNum}.frameNum;
        simSignal.timeSolution{jobNum}.sliceNum = spinModel.data{jobNum}.sliceNum;
        simSignal.timeSolution{jobNum}.numParts = spinModel.data{jobNum}.numParts;
        simSignal.timeSolution{jobNum}.partNum  = spinModel.data{jobNum}.partNum;
        
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
                '%s : done simulation %d/%d on GPU %d ',...
                functionName, jobNum, spinModel.numJobs,...
                spinModel.data{jobNum}.GPUindex);
            fprintf(fid, '  -- Elapsed Time   %.3fs', toc(tSimulation(jobId)));
            fprintf(fid, '\n');
        end
        
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