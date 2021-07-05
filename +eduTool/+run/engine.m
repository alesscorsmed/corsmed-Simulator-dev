function [simSignal] = engine( spinModel, pulseSequence, expControl)
%
% EDUTOOL.RUN.ENGINE
%
%	Runs the simulator.
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.engine';
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
expControl.estimatedRunTime = seconds(expControl.estimatedRunTime);
expControl.estimatedRunTime.Format = 'mm:ss';
eduToolMessage = sprintf('Pedal to the Metal... (Estimated Run Time %s sec)',...
    strrep(sprintf('%s',expControl.estimatedRunTime),':', ' min '));

% intialize progress to 5%
expControl.progress = 5;
eduTool.frontend.updateScannerStatus(expControl,eduToolMessage);
eduTool.frontend.updateExperimentProgress(expControl);

rxTime = pulseSequence.time(pulseSequence.rxSignal(:) > 0);
%% loop on the number of Jobs
for jobNum = 1:spinModel.numJobs
    
    %% report start of simulation
    if expControl.debug.debugMode
        tSimulation = tic();
        fprintf(fid, ...
            '\n\n%s : (%s) starting simulation %d/%d',...
            functionName, expControl.simulation.simulationEngine,...
            jobNum, spinModel.numJobs);
        fprintf(fid, '\n');
    end
    
    %% extract the model from the job queue
    model = spinModel.data{jobNum}.model;
    GPUindex = spinModel.data{jobNum}.GPUindex;
    %% run the simulation
    timeSolution = [];
    [timeSolution] = simulator.runSimulation( ...
        model, pulseSequence, expControl, timeSolution, GPUindex );
    
    %% store the solution
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
    simSignal.timeSolution{jobNum}.partNum = spinModel.data{jobNum}.partNum;
    
    %% update progress bar (simulation assumed 90%)
    expControl.progress = expControl.progress +...
        90/spinModel.numJobs;
   eduTool.frontend.updateScannerStatus(expControl);
    %% check status
    eduTool.frontend.checkExperimentStatus(...
        expControl,expControl.experimentID);
    
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