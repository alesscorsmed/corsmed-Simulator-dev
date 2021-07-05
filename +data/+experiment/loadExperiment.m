function [expControl,acquisition,inputSource] = loadExperiment(...
    inputSource,application,approach,mode)
%
% DATA.EXPERIMENT.LOADEXPERIMENT
%
%	Function that loads expControl and acquisition structs 
%   with data from an App-dependent source
%
% INPUT
%   inputSource   can be different types depending on the application
%   application   string with application type
%
% OUTPUT
%   expControl   expControl structure
%   acquisition  acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%

functionName = 'data.experiment.loadExperiment';
if (nargin < 2)
    ME = MException('Load:invalidArgCount',...
        '%s : Wrong argument count',functionName);
    throw(ME);
elseif nargin<3
    approach    = '';
    mode        = '';
end

%% initialize the structure with defaults
[expControl]    = data.experiment.initializeExpControl(approach);
[acquisition]   = data.experiment.initializeAcquisition();

%% check the app and load the data
% try
    switch lower(application)

        case lower('edutool')

            switch lower(approach)

                case lower('jsonstandalone')
                    
                    % update the app, the approach and the simulator mode
                    expControl.application  = 'edutool';
                    expControl.approach     = 'jsonstandalone';
                    expControl.mode         = mode;
                    
                    %% edutool-jsonstandalone
                    % get from redis the json file that holds all the experiments
                    [experimentsStr,~,~]    = tools.redis.redisGet(...
                        inputSource.redis.R,...
                        inputSource.redis.keys.experimentsRedisKey);
                    
                    %strsplit(strtrim(experimentsStr{1,1}),',');

                    if isempty(strtrim(experimentsStr{1,1}))
                        expControl  = '';
                        acquisition = '';
                        return;
                    end
                    
                    experimentArray         = jsondecode(experimentsStr{1,1}); 
                    
                    % the simulator continues only if 'kernel' or 'confirmed'
                    % experiments are available in the 'EXPERIMENT_{InstanceID}'
                    % redis key
                    idx         = find(ismember({experimentArray.status},'kernel') |...
                        ismember({experimentArray.status},'confirmed'));
                    
                    % the simulator should not continue if an experiment 
                    % with "cancelled-error" or "confirm" status is pending
                    idpending   = find(ismember({experimentArray.status},'cancelled-error') |...
                        ismember({experimentArray.status},'confirm'));
                    
                    if isempty(idx) || (~isempty(idpending) && ...
                            (str2double(experimentArray(idpending(1,1)).experiment_id) <  ...
                            str2double(experimentArray(idx(1,1)).experiment_id)))
                        expControl  = '';
                        acquisition = '';
                        return;
                    else
                        nextExpID = ...
                            str2double(experimentArray(idx(1,1)).experiment_id);
                        
                        if strcmp(experimentArray(idx(1,1)).status,'confirmed')
                            inputSource.latestConfirmedExperimentID = nextExpID;
                        else
                            inputSource.latestConfirmedExperimentID = 0;
                        end
                        
                        inputSource.latestExperimentID              = nextExpID;
                        inputSource.redis.experimentArray           = experimentArray;
                    end

                    % Get the redis key with all the data of the experiment
                    experimentRedisKey          = ...
                        [inputSource.redis.keys.experimentsRedisKey,...
                        '_',num2str(nextExpID)];
                    experimentUpdatesRedisKey   = ...
                        [inputSource.redis.keys.experimentsRedisKey,...
                        '_',num2str(nextExpID),'_UPDATES'];
                    experimentInfoRedisKey   = ...
                        [inputSource.redis.keys.experimentsRedisKey,...
                        '_',num2str(nextExpID),'_INFO'];
                    experimentResultsRedisKey   = ...
                        [inputSource.redis.keys.experimentsRedisKey,...
                        '_',num2str(nextExpID),'_RESULTS'];

                    [experimentData,inputSource.redis.R,redisReport] = ...
                        tools.redis.redisGetJsonWrapper(inputSource.redis.R,...
                        experimentRedisKey);
                    
                    % update timestamp & filetimestamp to be equal to the 
                    % one brought by jsonRedis
                    erasedChars = ["-",":","T"];
                    expControl.fileTimeStamp    = erase(extractBefore(...
                        experimentData.experiment.updatedAt,'.'),erasedChars);
                    expControl.timeStamp        = strrep(extractBefore(...
                        experimentData.experiment.updatedAt,'.'),'T',' ');

                    % If there is no key, or the key is empty, report an error for 
                    % this experiment and return
                    if isempty(experimentData) && ...
                            (strcmp(redisReport,'ERROR - NONEXISTANT KEY') || ...
                            strcmp(redisReport,'EMPTY KEY'))

                        errorMessage = fprintf(1,...
                            '%s: Error with the redis key for experiment %d\n',...
                            functionName,nextExpID);
                        
                        expControl  = '';
                        acquisition = '';

                        edutool.frontend.updateExperimentProgress(inputSource,...
                            '','error',errorMessage)
                        
                        ME = MException('eduTool:wrongArgCount',...
                            '%s : Error with the redis key for experiment %d',...
                            functionName,nextExpID);
                        throw(ME);
                    end

                    inputSource.redis.keys.experimentRedisKey           = ...
                        experimentRedisKey;
                    inputSource.redis.keys.experimentUpdatesRedisKey    = ...
                        experimentUpdatesRedisKey;
                    inputSource.redis.keys.experimentResultsRedisKey    = ...
                        experimentResultsRedisKey;
                    inputSource.redis.keys.experimentInfoRedisKey       = ...
                        experimentInfoRedisKey;

                    experimentData.redis =  inputSource.redis;

                     % for the experiment Control
                    [expControl] = services.setup.loadExpControl(...
                        expControl,experimentData,inputSource);

                    % for the acquisition
                    [acquisition] = services.setup.loadAcquisition(...
                        acquisition,experimentData);

                    % tissue Properties
                    [tissueData] = services.setup.loadTissueData(...
                        experimentData);

                    % upgraded expControl for return of tissueData
                    expControl.tissueData = tissueData;

                    % update expControl with the latest experiment ID
                    expControl.latestExperimentID = ...
                        inputSource.latestExperimentID;
                    
                    if ~isempty(acquisition) && isfield(acquisition.data,'generateB0map') 
                        expControl.b0map = acquisition.data.generateB0map;
                    end

                case lower('distributed')
                    
                    % update the app, the approach and the simulator mode
                    expControl.application  = 'edutool';
                    expControl.approach     = 'distributed';
                    expControl.mode         = '';
                    
                    %% edutool-distributed
                    % json file generated from the main API
                    % need to parse the json file to convert data

    %                 load('20210125_inputs_efsKernels.mat')
                    experimentData = inputSource;                

                    % for the experiment Control
                    [expControl] = services.setup.loadExpControl(...
                        expControl,experimentData);

                    % for the acquisition
                    [acquisition] = services.setup.loadAcquisition(...
                        acquisition,experimentData);

                    % tissue Properties
                    [tissueData] = services.setup.loadTissueData(...
                        experimentData);

                    % upgraded expControl for return of tissueData
                    expControl.tissueData = tissueData;

                case lower('expJson')
                    
                    % update the app, the approach and the simulator mode
                    expControl.application  = 'edutool';
                    expControl.approach     = 'expJson';
                    expControl.mode         = '';
                    
                    %% edutool-expjson
                    % json file with experiment data 

                    % load data from a valid json file inputSource
                    fid = fopen(inputSource,'r');
                    experimentData = jsondecode(fread(fid,inf,'*char').');
                    fclose(fid);

                    % assign data by deep copy of experiment structures
                    % this will avoid crashes in case there are 
                    % new fields that are not in the experiment data
                    expControl  = tools.misc.deepCopyStruct( ...
                        expControl, experimentData.expControl);
                    acquisition = tools.misc.deepCopyStruct( ...
                        acquisition, experimentData.acquisition);

                    % change the app
                    expControl.application = 'edutool';

                    % update folder system
                    % update kernel path
                    expControl.simulation.kernelFolder = ...
                        sprintf('%s/PROJECTS/%s/kernels/',...
                        expControl.folderSystem.baseFolder, expControl.application);
                    % default kernel
                    expControl.simulation.kernelPtx = sprintf('%s%s.ptx', ...
                        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
                    % upgrade with more info
                    expControl.folderSystem.errorFolder = ...
                        sprintf('%s/ERRORS/%s/', ...
                        expControl.folderSystem.baseFolder, expControl.application);
                    expControl.folderSystem.statsFolder = ...
                        sprintf('%s/STATS/%s/', ...
                        expControl.folderSystem.baseFolder, expControl.application);
                    expControl.folderSystem.experimentFolder = ...
                        sprintf('%s/RESULTS/%s/user_%d/', ...
                        expControl.folderSystem.baseFolder, ...
                        expControl.application, expControl.userID);
                    if ~isfolder(expControl.folderSystem.experimentFolder)
                        system(sprintf('mkdir %s',expControl.folderSystem.experimentFolder));
                    end

                otherwise
                    %% 'eduTool', defined by sessionData with DB queries   
                    
                    if isfield(inputSource,'mode') && ...
                            strcmp(inputSource.mode,'test')
                        
                        % load expControl and acquisition from external
                        % file
                        if isempty(inputSource.jsontestfile)
                            load('inputs/edutool-standalone-test_20210211_A6.mat')
                            
                            % update the app, the approach and the simulator mode
                            expControl.application  = 'edutool';
                            expControl.approach     = 'standalone';
                            expControl.mode         = 'test';
                        else
                            fid = fopen(inputSource.jsontestfile,'r');
                            experimentData = jsondecode(fread(fid,inf,'*char').');
                            fclose(fid);

                            % assign data by deep copy of experiment structures
                            % this will avoid crashes in case there are 
                            % new fields that are not in the experiment data
                            expControl  = tools.misc.deepCopyStruct( ...
                                expControl, experimentData.expControl);
                            acquisition = tools.misc.deepCopyStruct( ...
                                acquisition, experimentData.acquisition);
                            
                            % update the user id
                            expControl.userID = inputSource.userID;
                            
                            % update the app, the approach and the simulator mode
                            expControl.application  = 'edutool';
                            expControl.approach     = 'standalone';
                            expControl.mode         = 'test';

                            % update folder system
                            % update kernel path
                            expControl.simulation.kernelFolder = ...
                                sprintf('%s/PROJECTS/%s/kernels/',...
                                expControl.folderSystem.baseFolder, expControl.application);
                            % default kernel
                            expControl.simulation.kernelPtx = sprintf('%s%s.ptx', ...
                                expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
                            % upgrade with more info
                            expControl.folderSystem.errorFolder = ...
                                sprintf('%s/ERRORS/%s/', ...
                                expControl.folderSystem.baseFolder, expControl.application);
                            expControl.folderSystem.statsFolder = ...
                                sprintf('%s/STATS/%s/', ...
                                expControl.folderSystem.baseFolder, expControl.application);
                            expControl.folderSystem.experimentFolder = ...
                                sprintf('%s/RESULTS/%s/user_%d/', ...
                                expControl.folderSystem.baseFolder, ...
                                expControl.application, expControl.userID);
                            if ~isfolder(expControl.folderSystem.experimentFolder)
                                system(sprintf('mkdir %s',expControl.folderSystem.experimentFolder));
                            end
                        end
                        
                        % update the connection to db
                        expControl.connLocalDB  = inputSource.connLocalDB;
                        
                        % add a random pause for stress-tests on multiple
                        % tests
                        pauseTime = randsample(0:inputSource.testPauseTime,1);
                        pause(pauseTime)
                                                
                    else
                        % inputSource should be a sessionData with DB conn
                        [expControl] = eduTool.setup.dbQueryExpControl(...
                            expControl,inputSource);

                        % for the acquisition, we need expControl with DB conn as source
                        [acquisition] = eduTool.setup.dbQueryAcquisition(...
                            acquisition,expControl);
            
                        if ~isempty(acquisition) && isfield(acquisition.data,'generateB0map') 
                            expControl.b0map = acquisition.data.generateB0map;
                        end
                    end
                    
            end

        case lower('sbr')

            %% specific for the SBR, with a pre-defined structure
            % input source is a sessionData struct
            sessionData = inputSource; 

            % fill the relevant fields
            expControl.application                  = application;
            expControl.debug.errorFolder            = sessionData.errorFolder;
            expControl.simulation.kernelFolder      = sessionData.kernelFolder;
            expControl.simulation.kernelPtx         = sessionData.kernelPtx;
            expControl.simulation.simulationEngine  = sessionData.simulationEngine;
            expControl.simulation.numGPUs           = sessionData.numGPU;
            expControl.simulation.threads           = sessionData.threads;
            expControl.simulation.blocks            = sessionData.blocks;
            expControl.debug.debugMode              = sessionData.debugMode;
            expControl.debug.debugFile              = sessionData.debugFile;

            % Experiment name
            expControl.name = sprintf('%s_%s', ...
                application, expControl.fileTimeStamp);

        otherwise

            %% 'eduTool', defined by sessionData with DB queries
            % inputSource should be a sessionData with DB conn
            [expControl] = eduTool.setup.dbQueryExpControl(...
                expControl,inputSource);

            % for the acquisition, we need expControl with DB conn as source
            [acquisition] = eduTool.setup.dbQueryAcquisition(...
                acquisition,expControl);
            
            if ~isempty(acquisition) && isfield(acquisition.data,'generateB0map') 
                expControl.b0map = acquisition.data.generateB0map;
            end
    end
% catch ME
%     fprintf(1,['Error in function %s() at line %d.',...
%         '\n Error Message: %s\n'], ....
%         ME.stack(1).name,ME.stack(1).line,...
%         ME.message);   
% end

%% if there is no DB connection: add empty for consistency
% avoids crashing when updating state/progress
if ~isfield(expControl,'connLocalDB') && ~strcmp(approach,'jsonstandalone')
    expControl.connLocalDB = [];
end

