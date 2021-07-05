function [batchResult] = runBatch(batchFolder,batchString,updateJason)
%
% EDUTOOL.RUNBATCH
%
%	Runs a batch of experiments that match the string
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.batch.runBatch';
%
if (nargin < 2)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end
if (nargin < 3) || isempty(updateJason)
    updateJason = 0;
end

%% initialize a pool with number of GPUs available
tParpool = tic();
edtPool = [];
gpuPool = [];
numGPUs = gpuDeviceCount;
if numGPUs > 0
    % initialize parPool if it does not exist
    if isempty(gcp('nocreate'))
        edtPool = parpool(numGPUs, 'IdleTimeout', Inf);
    else
        edtPool = gcp();
        % if wrong size, restart
        if edtPool.NumWorkers > numGPUs
            delete(edtPool);
            edtPool = parpool(numGPUs, 'IdleTimeout', Inf);
        end
    end
    % create a gpuPool with the GPU devices
    spmd
        % each worker selects its own gpu
        gpuPool = gpuDevice(labindex);
    end
else
    MException('EduTool:BadInstance', 'Instance has not available GPUs');
end
fprintf(1, '\n Parallel Pool initialized with %d Workers', edtPool.NumWorkers);
fprintf(1, '\n  Initialization Time  %.3fs', toc(tParpool));
fprintf(1, '\n');  


%% get the experiments to run
experimentList = dir(sprintf('%s/experiments/*.json',batchFolder));

%% find the files with any of the strings (+ as separator)
batchString = split(batchString,'+');
idxToRun = find(contains(lower({experimentList(:).name}),lower(batchString(:))));
courseID = -1;

%% allocate for the results
batchResult.numExperiments  = length(idxToRun);
batchResult.numPass         = 0;
batchResult.numFail         = 0;
% allocate for each experiment
batchResult.experiment{length(idxToRun)} = [];


%% loop on the files and simulate
for ii = 1:length(idxToRun)
    
    %% experiment to run
    experimentName = experimentList(idxToRun(ii)).name;
    experimentFile = sprintf('%s/experiments/%s',batchFolder,experimentName);
    
    %% add to Result
    batchResult.experiment{ii}.name         = experimentName;
    batchResult.experiment{ii}.jsonUpdate   = 0;
    
    try
        
        %% load experiment data from json file
        [expControl,acquisition] = data.experiment.loadExperiment(...
            experimentFile,'expJson');
        
        %% update acquisition with the computed Noise levels
        [acquisition] = noise.calculateNoiseLevel( acquisition, ...
            expControl.mrSystem.b0 );  
    
        %% If need to update the json file:
        %  re-write into json file (overwrite with new structs)
        if updateJason
            
            batchResult.experiment{ii}.jsonUpdate = 1;
            
            experimentData.acquisition = acquisition;
            experimentData.expControl  = expControl;
            experimentData = jsonencode(experimentData);
            fprintf(1, '\n');
            fprintf(1, '\n Saving experiment ... ');
            tJson = tic();
            % save experimentData to json file
            fid = fopen(experimentFile,'w');
            fwrite(fid, experimentData, 'char');
            fclose(fid);
            fprintf(1, ' done -- elapsed time %.3fs', toc(tJson));
            fprintf(1, '\n  File : %s', experimentFile);
            fprintf(1, '\n');
        end
        
        %% load anatomical model
        if (courseID ~= expControl.courseID) || isempty(anatomicalModel)
            [anatomicalModel] = data.models.initializeAnatomical(...
                expControl.courseID, ...
                expControl.folderSystem.anatomicalModelFolder, ...
                expControl.application);
        end
        
        %% load coil system
        if (courseID ~= expControl.courseID) || isempty(coilSystem)
            [coilSystem] = data.models.initializeCoils(...
                expControl.courseID, ...
                expControl.folderSystem.coilModelFolder);
            
            %% precompute SAR for tx coils
            [coilSystem] = coils.precomputeMRSafety(coilSystem, anatomicalModel);
            
        end
        
        % assign current course ID
        courseID  = expControl.courseID;
        
        %% Update Coil Selection
        [coilSystem] = coils.updateCoilSelection(...
            expControl, anatomicalModel, coilSystem );
        
        %% upgrade expControl with parPool and gpuPool
        [expControl,edtPool,gpuPool,numGPUs] = ...
            eduTool.multiGPU.setParallelization( ...
            expControl,edtPool,gpuPool,numGPUs);
        
        %% Start the Experiment
        tExperiment = tic();
        fprintf(1, '\n STARTING IMAGING EXPERIMENT %s', experimentName);
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
        fprintf(1, '\n IMAGING EXPERIMENT %s DONE', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);
        fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
        
        %% noted as pass
        batchResult.numPass               = batchResult.numPass + 1;
        batchResult.experiment{ii}.result = 'PASS';
        batchResult.experiment{ii}.note   = sprintf('Elapsed time : %.3fs', tExperiment);
        
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
        fprintf(1, '\n IMAGING EXPERIMENT %s FAILED', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);
        fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
   
        %% add as failed to the result
        batchResult.numFail               = batchResult.numFail + 1;
        batchResult.experiment{ii}.result = 'FAIL';
        batchResult.experiment{ii}.note   = ...
            sprintf('Error in function %s() at line %d. Message: %s',...
            ME.stack(1).name,ME.stack(1).line,ME.message);
        
    end

end
