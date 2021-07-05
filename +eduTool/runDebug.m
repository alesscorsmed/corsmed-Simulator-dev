%
% EDUTOOL.RUNDEBUG
%
%	Runs a experiment for debuging
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.runDebug';

%% define the experiment json file
% experimentFile = '/efs-mount-point/S20/ERRORS/edutool/20210125130448_U933_C6_E339_P360_UI255_ERROR.json';
% experimentFile = '/efs-mount-point/S20/ERRORS/edutool/20210126135613_U1139_C6_E276_P288_UI228_ERROR.json';
%experimentFile = '/efs-mount-point/S20/ERRORS/edutool/20210212094637_U1146_C7_E79_P80_UI76_ERROR.json';
experimentFile = '/efs-mount-point/S20/RESULTS/edutool/user_790/20210330103635_U790_C4_E338_P348_UI322_ERROR.json';

%% initialize a pool with number of GPUs available
tParpool = tic();
CPUperGPU = 2;
edtPool = [];
gpuPool = [];
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
fprintf(1, '\n Parallel Pool initialized with %d Workers', edtPool.NumWorkers);
fprintf(1, '\n  Initialization Time  %.3fs', toc(tParpool));
fprintf(1, '\n');

tExperiment = tic();
try
    
    %% load experiment data from json file
    [expControl,acquisition] = data.experiment.loadExperiment(...
        experimentFile,'edutool','expJson','debug');
    
    
    %% update acquisition with the computed Noise levels
    [acquisition] = noise.calculateNoiseLevel( acquisition, ...
        expControl.mrSystem.b0 );
    
    %% load anatomical model
    if ~exist('anatomicalModel', 'var') || isempty(anatomicalModel)
        [anatomicalModel] = data.models.initializeAnatomical(...
            expControl.courseID, ...
            expControl.folderSystem.anatomicalModelFolder, ...
            expControl.application);
    end
    
    %% load coil system
    if ~exist('coilSystem', 'var') ||isempty(coilSystem)
        [coilSystem] = data.models.initializeCoils(...
            expControl.courseID, ...
            expControl.folderSystem.coilModelFolder);
        
        %% precompute SAR for tx coils
        [coilSystem] = coils.precomputeMRSafety(coilSystem, anatomicalModel);
    end
    
    
    %% Update Coil Selection
    [coilSystem] = coils.updateCoilSelection(...
        expControl, anatomicalModel, coilSystem );
    
    %% update application and user
    expControl.application      = 'debug';
    expControl.userID           = '000';
    expControl.useOldSequence   = 0; % old pulse sequence does not work
    
    
    %% upgrade expControl with parPool and gpuPool
    [expControl,edtPool,gpuPool,numGPUs] = ...
        eduTool.multiGPU.setParallelization( ...
        expControl,edtPool,gpuPool,numGPUs);
    
    %% Start the Experiment
    tExperiment = tic();
    fprintf(1, '\n STARTING IMAGING EXPERIMENT %s', expControl.name);
    fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
    fprintf(1, '\n  User       ID : %d', expControl.userID);
    fprintf(1, '\n  Course     ID : %d', expControl.courseID);
    fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
    fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
    fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);
    fprintf(1, '\n');
    
    %% call the execution
    info4user = eduTool.run.experimentExecution(...
        anatomicalModel, coilSystem, acquisition, expControl);
    
    %% report experiment done
    tExperiment = toc(tExperiment);
    
    %% print result
    fprintf(1, '\n');
    fprintf(1, '\n IMAGING EXPERIMENT %s DONE', expControl.name);
    fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
    fprintf(1, '\n  User       ID : %d', expControl.userID);
    fprintf(1, '\n  Course     ID : %d', expControl.courseID);
    fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
    fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
    fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);
    fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
    fprintf(1, '\n');
    
    %% noted as pass
    pass = 1;
    
catch ME
    
    %% catch the error
    ME.identifier;
    ME.message;
    
    %% send error in connection to DB for backend
    errorMessage = sprintf(['%s - User (%d) - Exper (%d)',...
        '\n Error in function %s() at line %d.',...
        '\n Error Message: %s',...
        '\n Data for error replication saved in %s'], ...
        expControl.timeStamp, ...
        expControl.userID,...
        expControl.experimentID,...
        ME.stack(1).name,ME.stack(1).line,ME.message,...
        experimentFile);
    
    %% display error in cmd line
    tExperiment = toc(tExperiment);
    fprintf(1, '\n');
    fprintf(1, '\n ERROR: %s', errorMessage);
    fprintf(1, '\n');
    fprintf(1, '\n IMAGING EXPERIMENT %s FAILED', expControl.name);
    fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
    fprintf(1, '\n  User       ID : %d', expControl.userID);
    fprintf(1, '\n  Course     ID : %d', expControl.courseID);
    fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
    fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
    fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);
    fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
    fprintf(1, '\n');
    
    %% add as failed to the result
    pass = 0;
    
end
