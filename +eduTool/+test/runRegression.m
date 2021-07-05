function [testResult] = runRegression(regressionFolder,testString, testMode)
%
% EDUTOOL.TEST.RUNREGRESSION
%
%	Runs a batch of tests of experiments that match the string from the
%	subfolder /inputs in the regression folder.
%   If testMode = 'generation'
%       Generates the reference (golden) data and stores in the
%       subfolder /golden in the regression folder
%   If testMode = 'test'
%       Runs the experiment and compares the results with glodenr
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool:test:runRegression';
%
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end
if (nargin < 2)
    testString = []; % runs all available cases
end
if (nargin < 3) || isempty(testMode)
    testMode = 'test'; % test / generation
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
experimentList = dir(sprintf('%s/inputs/*.json',regressionFolder));

%% find the files with any of the strings (+ as separator)
if isempty(testString) || strcmpi(testString, 'all')
    idxToRun = 1:numel(experimentList);
else
    testString = split(testString,'+');
    idxToRun = find(contains(lower({experimentList(:).name}),lower(testString(:))));
end
courseID = -1;

%% allocate for the results
testResult.numExperiments  = length(idxToRun);
testResult.numPass         = 0;
testResult.numFail         = 0;
% allocate for each experiment
testResult.experiment{length(idxToRun)} = [];


%% loop on the files and simulate
for ii = 1:length(idxToRun)
    
    %% experiment to run
    [~,experimentName,ext] = fileparts(experimentList(idxToRun(ii)).name);
    experimentFile = sprintf('%s/inputs/%s%s',regressionFolder,experimentName,ext);
    
    %% add to Result
    testResult.experiment{ii}.name	= experimentName;
    testResult.experiment{ii}.file  = experimentFile;
    
    timeStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
    % change time stamp into file format yyyymmddHHMMSS
    timeStamp = strrep(timeStamp,'-','');
    timeStamp = strrep(timeStamp,':','');
    timeStamp = strrep(timeStamp,' ','');
    
    %% prepare the testControl for the experiment
    testControl = [];
    testControl.regFolder   = regressionFolder; % regression folder
    testControl.testName    = experimentName;   % name of file
    if strcmpi(testMode, 'generation')
        testControl.testGen     = 1; % generation of data
        testControl.testRun     = 0; % run the regression test
    else
        testControl.testGen     = 0; % generation of data
        testControl.testRun     = 1; % run the regression test
    end
    testControl.status.pass = 1; % flag
    testControl.version     = 'v2.9.19';
    testControl.timeStamp   = timeStamp;
    
    try
        
        %% load experiment data from json file
        [expControl,acquisition] = eduTool.test.loadExperimentJson(...
            experimentFile);
                
        %% update acquisition with the computed Noise levels
        [acquisition] = noise.calculateNoiseLevel( acquisition, ...
            expControl.mrSystem.b0 );  
        
        %% load anatomical model
        if (courseID ~= expControl.courseID) || isempty(anatomicalModel)
            [anatomicalModel] = data.models.initializeAnatomical(...
                expControl.anatomicalID, ...
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
        % remove any debugging user
        expControl.application      = 'regression';
        expControl.approach         = 'regression-testing';
        expControl.userID           = '000';
        expControl.debug.debugMode  = 0;
        % update the results folder
        expControl.folderSystem.experimentFolder = ...
            sprintf('%s/result/',regressionFolder);
        
        %% Update Coil Selection
        [coilSystem] = coils.updateCoilSelection(...
            expControl, anatomicalModel, coilSystem );
        
        %% upgrade expControl with parPool and gpuPool
        [expControl,edtPool,gpuPool,numGPUs] = ...
            eduTool.multiGPU.setParallelization( ...
            expControl,edtPool,gpuPool,numGPUs);
        
        %% Start the Experiment
        tExperiment = tic();
        fprintf(1, '\n STARTING: %s ... \n\n', experimentName);
        
        %% call the execution
        [~, testControl] = eduTool.run.experimentExecution(...
            anatomicalModel, coilSystem, acquisition, ...
            expControl, [], testControl);      
        
        %% report experiment done
        tExperiment = toc(tExperiment);
        
        %% print result
        fprintf(1, '\n\n DONE :  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
        if testControl.status.pass
            %% noted as pass
            testResult.numPass               = testResult.numPass + 1;
            testResult.experiment{ii}.result = 'PASS';
        else
            %% noted as fail
            testResult.numFail               = testResult.numFail + 1;
            testResult.experiment{ii}.result = 'FAIL';
        end
        testResult.experiment{ii}.note   = sprintf('Elapsed time : %.3fs', tExperiment);
        % print details
        statusFields = fieldnames(testControl.status);
        % loop on the field names and call recursive for structs, otherwise assign
        for ff = 1:length(statusFields)
            testResult.experiment{ii}.report.(statusFields{ff}) = ...
                testControl.status.(statusFields{ff});
        end
        
    catch ME
        
        %% catch the error
        ME.identifier;
        ME.message;

        %% send error in connection to DB for backend
        errorMessage = sprintf(['CRASH: ',...
            'function %s() at line %d. Message: %s'], ...
            ME.stack(1).name, ME.stack(1).line, ME.message);

        %% display error in cmd line
        fprintf(1, '\n FAIL : %s', errorMessage);
        fprintf(1, '\n');
   
        %% add as failed to the result
        testResult.numFail               = testResult.numFail + 1;
        testResult.experiment{ii}.result = 'CRASH';
        testResult.experiment{ii}.note   = ...
            sprintf('Error in function %s() at line %d. Message: %s',...
            ME.stack(1).name,ME.stack(1).line,ME.message);
        
    end

end
