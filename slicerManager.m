function slicerManager(redisAddress,redisPort,...
    redisRequestKey,redisResponseKey,redisUpdatesKey,calcIDkey,...
    anatModelsIDstring,anatomicalFolder,coilsFolder)
%
% SERVICES.SLICERMANAGER
%
%	Mock up of the SLICER manager service
%
% INPUT
%
% redisAddress (string)
% redisPort (string)
% redisRequestKey (string): the redis key that holds the name of the key
%   that keeps the experiment data in json format
% redisResponseKey (string): the redis key where this service will store 
%   the name of the key thats keep all the results from the calculator
% redisUpdatesKey (string): the redis key where this service will store 
%   status of the experiment and its progress
% calcIDkey: the redis key that holds the unique name of the variables 
%   stored in redis for future use by the simulator and the reconstructor
% anatModelsIDstring (string): it holds all the IDs of the available
%   anatomical models. Example: '[4,6,7]' or '[1:7]'
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'slicerManager';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% Establish connection with Redis
% Connect to a redis server
redis.address   = redisAddress;
redis.port      = redisPort;

R = tools.redis.redisEstablishConnection(redis.address,redis.port);

%% Load all models from local filesystem
timerFilesLoad  = tic;
anatIDs         = eval(anatModelsIDstring);

% Preallocate memory
maxAnatID                       = max(anatIDs);
anatomicalModelAll{1,maxAnatID} = [];
coilSystemAll{1,maxAnatID}      = [];
for anatID = anatIDs
    tTotal = tic;
    %% load anatomical model
    fprintf(1,'Loading the anatomical model %i\n',anatID);
    anatomicalModelAll{1,anatID} = data.models.initializeAnatomical(anatID,...
        [anatomicalFolder,filesep],'edutool');
    
    %% load coil system
    fprintf(1,'Coil calibration and setup for anatomical model %i\n',anatID);
    coilSystemAll{1,anatID} = data.models.initializeCoils(anatID,...
        coilsFolder);
    
    %% precompute SAR for tx coils
    fprintf(1,'MR safety checks for anatomical model %i\n',anatID);
    coilSystemAll{1,anatID} = coils.precomputeMRSafety(coilSystemAll{1,anatID},...
        anatomicalModelAll{1,anatID});
    
    %% report
    tTotal = toc(tTotal);
    fprintf(1, '\n%s : done for course ID %d',...
        functionName, anatID );
    fprintf(1, '\n  Scanner ready' );
    fprintf(1, '\n    Active coil     %s', coilSystemAll{1,anatID}.activeRx);
    fprintf(1, '\n    Patient name    %s', anatomicalModelAll{1,anatID}.name);
    fprintf(1, '\n  Booting Time      %.3fs', tTotal);
    fprintf(1, '\n\n');
end
timerFilesLoadTotal = toc(timerFilesLoad);
fprintf(1, '\n  TOTAL BOOTING TIME      %.3fs', timerFilesLoadTotal);

%% While loop that checks on the redisRequestKey
while 1
    
    try
        % Check if the key exists and it is not empty
        [jsonRedis,~,redisReport] = tools.redis.redisGet(R,redisRequestKey);

        % In case MATLAB has disconnected from redis, try to reconnect
        if strcmp(redisReport,'ERROR - NO CONNECTION')
            R = tools.redis.redisEstablishConnection(redis.address,redis.port);
        end

        % If there is no key, or the key is empty, continue
        if isempty(jsonRedis) || strcmp(redisReport,'ERROR - NONEXISTANT KEY')
            continue;
        elseif iscell(jsonRedis)
            if isempty(jsonRedis{1,1})
                continue;
            end        
        end
        
        % Fetch the uniqueID of this experiment
        [uniqueIDcell,~,~] = tools.redis.redisGet(R,calcIDkey);
        uniqueID = uniqueIDcell{1,1};

        %% Fetch data from Redis
        tRedis      = tic();
        
        % Extract information for running the experiment
        experimentData = jsondecode(jsonRedis{1,1});
        
        % Ask MOP to add the following fields in json @@@
        experimentData.application = 'apijson';
        
        fprintf(1, '\n  Redis retrieval time  : %.3fs', toc(tRedis));
        
        
        %% Initialize experiment
        tTotal = tic();
        fprintf(1, '\n STARTING SLICE MANAGER for redis-key: %s',redisRequestKey);
        fprintf(1, '\n');
        
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'update','Experiment is running','1');
        % Returns the necessary information for the experiment from FrontEnd
        [~, expControl, acquisition] = ...
            eduTool.setup.fetchExperimentData(experimentData);
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'update','Experiment is running','3');

        expControl.redis    = redis;
        
        %% Fields that we have to update
        expControl.versionNum                           = 'v20210120_S20';

        %% TEMP - ask MOP to add these fields in the json string stored in Redis
        expControl.anatomicalID                         = 6;
        expControl.application                          = 'edutool';
        expControl.useOldSequence                       = 0;
        expControl.applicationType                      = experimentData.application;
        %% TEMP - end @@@
        
        anatomicalModel = anatomicalModelAll{1,expControl.anatomicalID}; 
        coilSystem      = coilSystemAll{1,expControl.anatomicalID};        

        %% Update Coil Selection
        [coilSystem] = coils.updateCoilSelection(...
            expControl, anatomicalModel, coilSystem );
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'update','Experiment is running','5');

        %% CALL THE SLICER EXECUTION
        [spinModel, pulseSequence, motionModel, sarReport, imageData, reconData, ...
            expControl] = eduTool.run.slicer( anatomicalModel, coilSystem, ...
            expControl.mrSystem, acquisition, expControl );
        
%         [simJobs, spinModel, pulseSequence, motionModel, sarReport, ...
%             acquisition, expControl] = eduTool.run.slicerV3(...
%             anatomicalModel, coilSystem, mrSystem, acquisition, expControl );
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'update','Experiment is running','8');

        %% SAVE DATA TO REDIS
        % Establish redis connection using the MATLAB built-in function
        c = tools.redisMatlab.redisEstablishConnectionMatlab(...
            redis.address,redis.port);
        
        % Add spinModel in redis and verify that it is available in redis
        spinModelKey = [uniqueID,'_spinModel'];
        put(c,spinModelKey,spinModel);
        if ~isKey(c,spinModelKey)
            ME = MException('eduTool:spinModelRedisSetFailure',...
                '%s : could not be stored in redis',spinModelKey);
            throw(ME);
        end
        
        % Add pulseSequence in redis and verify that it is available in redis
        pulseSequenceKey = [uniqueID,'_pulseSequence'];
        put(c,pulseSequenceKey,pulseSequence);
        if ~isKey(c,pulseSequenceKey)
            ME = MException('eduTool:pulseSequenceRedisSetFailure',...
                '%s : could not be stored in redis',pulseSequenceKey);
            throw(ME);
        end
        
        % Add motionModel in redis and verify that it is available in redis
        motionModelKey = [uniqueID,'_motionModel'];
        put(c,motionModelKey,motionModel);
        if ~isKey(c,motionModelKey)
            ME = MException('eduTool:motionModelRedisSetFailure',...
                '%s : could not be stored in redis',motionModelKey);
            throw(ME);
        end
        
        % Add acquisition in redis and verify that it is available in redis
        acquisitionKey = [uniqueID,'_acquisition'];
        put(c,acquisitionKey,acquisition);
        if ~isKey(c,acquisitionKey)
            ME = MException('eduTool:acquisitionRedisSetFailure',...
                '%s : could not be stored in redis',acquisitionKey);
            throw(ME);
        end
        
        % Add expControl in redis and verify that it is available in redis
        expControlKey = [uniqueID,'_expControl'];
        put(c,expControlKey,expControl);
        if ~isKey(c,expControlKey)
            ME = MException('eduTool:expControlRedisSetFailure',...
                '%s : could not be stored in redis',expControlKey);
            throw(ME);
        end
        
        % Add sarReport in redis and verify that it is available in redis
        sarReportKey = [uniqueID,'_sarReport'];
        put(c,sarReportKey,sarReport);
        if ~isKey(c,sarReportKey)
            ME = MException('eduTool:sarReportRedisSetFailure',...
                '%s : could not be stored in redis',sarReportKey);
            throw(ME);
        end
        
        % Add imageData in redis and verify that it is available in redis
        imageDataKey = [uniqueID,'_imageData'];
        put(c,imageDataKey,imageData);
        if ~isKey(c,imageDataKey)
            ME = MException('eduTool:sarReportRedisSetFailure',...
                '%s : could not be stored in redis',imageDataKey);
            throw(ME);
        end
        
        % Add reconData in redis and verify that it is available in redis
        reconDataKey = [uniqueID,'_reconData'];
        put(c,reconDataKey,reconData);
        if ~isKey(c,reconDataKey)
            ME = MException('eduTool:reconDataRedisSetFailure',...
                '%s : could not be stored in redis',reconDataKey);
            throw(ME);
        end
        
        % save the simulation jobs independently and store the jobs in the 
        % redisResponseKey 
        redisResponseKeyValue = '';
        for jobNum = 1:spinModel.numJobs
            % extract and save the data
            simJob = spinModel.data{jobNum};
            simJobKey = [uniqueID,'_simJob',num2str(jobNum)];
            put(c,simJobKey,simJob);
            redisResponseKeyValue = strcat(redisResponseKeyValue,',',simJobKey);
            if ~isKey(c,simJobKey)
                ME = MException('eduTool:simJobRedis',...
                    '%s : could not be stored in redis',simJobKey);
                throw(ME);
            end
        end
        
        %% Delete redisRequestKey from redis
        [~, statusDel] = tools.redis.redisDEL(R,redisRequestKey);
        if ~strcmp(statusDel,'OK')
            ME = MException('eduTool:delKeyRedisFailure',...
                '%s : could not be deleted from redis',redisRequestKey);
            throw(ME);
        end
        
        %% Fill in the redisResponseKey with the redisResponseKeyValue
        [~, statusSet] = tools.redis.redisSet(R,redisResponseKey,...
            redisResponseKeyValue(2:end));
        if ~strcmp(statusSet,'OK')
            ME = MException('eduTool:setRedisResponseKeyFailure',...
                '%s : could not be set in redis',redisResponseKey);
            throw(ME);
        end
        
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'update','Experiment is running','10');

        %% REPORT OK
        fprintf(1, '\n');
        fprintf(1, '\n SLICE MANAGER on %s DONE', redisRequestKey);
        fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
        fprintf(1, '\n');

        %%
    catch ME
        %% send error back to UI
        errorMessage = sprintf(['Error in function %s() at line %d.',...
            '\n Error Message: %s'], ....
            ME.stack(1).name,ME.stack(1).line,...
            ME.message);
        % errorMessage = tools.printErrorMessage(expControl,ME); 
        
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'error',errorMessage);               
        
        tools.redis.redisDisconnect(R);
        
        disp(errorMessage)
       
    end
end