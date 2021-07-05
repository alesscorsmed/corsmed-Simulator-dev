function [simSignal] = engineService(spinModel,pulseSequence,motionModel,...
    expControl)
%
% EDUTOOL.RUN.ENGINESERVICE
%
%	Runs the simulator.
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.engineService';
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
expControl.estimatedRunTime.Format = 'mm:ss';
eduToolMessage = sprintf('Pedal to the Metal... (Estimated Run Time %s sec)',...
    strrep(sprintf('%s',expControl.estimatedRunTime),':', ' min '));
eduTool.frontend.updateScannerStatus(expControl,eduToolMessage);

% intialize progress to 10%
expControl.progress = 10;
rxTime = pulseSequence.time(pulseSequence.rxSignal(:) > 0);

%% define jobs to run
if isfield(expControl,'applicationType') ...
    && strcmp(expControl.applicationType,'apijson')
        jobs = expControl.numJob;
else
    jobs = 1:spinModel.numJobs;
end

%% loop on the number of jobs
for jobNum = jobs
    
    %% report start of simulation
    if expControl.debug.debugMode
        tSimulation = tic();
        fprintf(fid, ...
            '\n\n%s : starting simulation %d/%d',...
            functionName, jobNum, spinModel.numJobs);
        fprintf(fid, '\n');
    end
    
    %% extract the model from the job queue
    model = spinModel.data{jobNum}.model;
    GPUindex = spinModel.data{jobNum}.GPUindex;
    %% run the simulation
    timeSolution = [];
    [timeSolution] = simulator.runSimulation( ...
        model, pulseSequence, motionModel, expControl, timeSolution, GPUindex );
    
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
    
    %% report end of experiment
    if expControl.debug.debugMode
        fprintf(fid, ...
            '\n%s : done simulation %d/%d',...
            functionName, jobNum, spinModel.numJobs);
        fprintf(fid, '\n  Elapsed Time   %.3fs', toc(tSimulation));
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
eduTool.frontend.updateScannerStatus(expControl);

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Number of Jobs    %d', spinModel.numJobs);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end