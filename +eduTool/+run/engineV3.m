function [simSignal] = engineV3( spinModel, pulseSequence, motionModel, ...
    expControl)
%
% EDUTOOL.RUN.ENGINE
%
%	Runs the simulator, using spinTwin module
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.engineV3';
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
expControl.estimatedRunTime = seconds(expControl.estimatedRunTime);
expControl.estimatedRunTime.Format = 'mm:ss';
eduToolMessage = sprintf('Pedal to the Metal... (Estimated Run Time %s sec)',...
    strrep(sprintf('%s',expControl.estimatedRunTime),':', ' min '));

% intialize progress to 5%
expControl.progress = 5;
eduTool.frontend.updateScannerStatus(expControl,eduToolMessage);
eduTool.frontend.updateExperimentProgress(expControl);

%% extract simControl and dbgControl data from expControl
dbgControl.mode = expControl.debug.debugMode;
dbgControl.file = expControl.debug.debugFile;
% simulation control
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

rxTime = pulseSequence.time(pulseSequence.rxSignal(:) > 0);
%% loop on the number of Jobs
for jobNum = 1:spinModel.numJobs
    
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
    %% run the simulation
    [timeSolution, ~] = spinTwin.runSimulation( ...
        model, pulseSequence, motionModel, simControl, dbgControl);
    
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