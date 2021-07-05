function slicerManager3(redisAddress,redisPort,...
    redisRequestKey,redisResponseKey,redisUpdatesKey,calcIDkey,...
    anatModelsIDstring,anatomicalFolder,coilsFolder,breakLoop,dummyData)
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
% if (nargin < 1)
%     ME = MException('eduTool:wrongArgCount',...
%         '%s : wrong argument count',functionName);
%     throw(ME);
% end
disp('Calculator just started')
if (nargin < 10)
    breakLoop = 0;
    dummyData = 0;
end

%% Establish connection with Redis
% Connect to a redis server
disp('Trying to connect to redis')
redis.address       = redisAddress;
redis.port          = redisPort;
inputBufferSize     = 10000000;
outputBufferSize    = 10000000;

R = tools.redis.redisEstablishConnection(redis.address,redis.port,...
    inputBufferSize,outputBufferSize);

%% Load all models from local filesystem
disp('Start loading models:')
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
disp('Going into the while loop')
while 1
    
    try
        % Check if the key exists and it is not empty
%         [jsonRedis,~,redisReport] = tools.redis.redisGet(R,redisRequestKey);
        
        [experimentData,R,redisReport] = ...
            tools.redis.redisGetJsonWrapper(R,redisRequestKey);

        % If there is no key, or the key is empty, continue
        if isempty(experimentData) && ...
                (strcmp(redisReport,'ERROR - NONEXISTANT KEY') || ...
                strcmp(redisReport,'EMPTY KEY'))
            continue;     
        end
        
        % Fetch the uniqueID of this experiment
        [uniqueIDcell,~,~] = tools.redis.redisGet(R,calcIDkey);
        uniqueID = strtrim(uniqueIDcell{1,1});

        %% Fetch data from Redis
        tRedis      = tic();
        
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
        if dummyData
            load('20210125_inputs_efsKernels.mat')
        else
            [~, expControl, acquisition] = ...
                eduTool.setup.fetchExperimentData(experimentData);
        end
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
        
        tools.updateJsonProgress(R,redisUpdatesKey,...
            'update','Experiment is running','8');

        %% SAVE DATA TO REDIS
        
        pulseSequenceKey    = strcat(uniqueID,'_pulseSequence');
        motionModelKey      = strcat(uniqueID,'_motionModel');
        expControlKey       = strcat(uniqueID,'_expControl');
        sarReportKey        = strcat(uniqueID,'_sarReport');
        imageDataKey        = strcat(uniqueID,'_imageData');
        reconDataKey        = strcat(uniqueID,'_reconData');
        totalJobsKey        = strcat(uniqueID,'_totalJobs');
        
        redisResponseKeyValue = '';
        for jobNum = 1:spinModel.numJobs
            simJobKey{jobNum} = [uniqueID,'_simJob',num2str(jobNum)];
            redisResponseKeyValue = strcat(redisResponseKeyValue,',',...
                simJobKey{jobNum});
        end
        
        tTransferRedis = tic;
        fprintf(1,'\n  Storing key (json): %s',pulseSequenceKey)
        R = tools.redis.redisSetJsonWrapper(R,pulseSequenceKey,pulseSequence);
        fprintf(1,'\n  Storing key (json): %s',motionModelKey)
        R = tools.redis.redisSetJsonWrapper(R,motionModelKey,motionModel);

        % TEMPORARY FIX @@@
        reconData.encoding.operators        = [];
        expControl.estimatedRunTime         = 0;
        expControl.simulation.kernelFolder  = '/efs-mount-point/S20/PROJECTS/edutool/kernels';
        expControl.simulation.kernelPtx     = '/efs-mount-point/S20/PROJECTS/edutool/kernels/v26_sm70.ptx';

        
        fprintf(1,'\n  Storing key (json): %s',expControlKey)
        R = tools.redis.redisSetJsonWrapper(R,expControlKey,expControl);        
        fprintf(1,'\n  Storing key (json): %s',sarReportKey)
        R = tools.redis.redisSetJsonWrapper(R,sarReportKey,sarReport);
        fprintf(1,'\n  Storing key (json): %s',imageDataKey)
        R = tools.redis.redisSetJsonWrapper(R,imageDataKey,imageData);
        fprintf(1,'\n  Storing key (json): %s',reconDataKey)
        R = tools.redis.redisSetJsonWrapper(R,reconDataKey,reconData);

        % Flush data from input and output buffer of R
        flushinput(R)
        flushoutput(R)
        for jobNum = 1:spinModel.numJobs
            simJob = spinModel.data{jobNum};
            fprintf(1,'\n  Storing key (json): %s',simJobKey{jobNum})
            R = tools.redis.redisSetJsonWrapper(R,simJobKey{jobNum},simJob);
        end
        
        fprintf(1,'\n  Storing key: %s',totalJobsKey)
        [~, statusSet] = tools.redis.redisSet(R,totalJobsKey,...
            num2str(spinModel.numJobs));

        tTransferRedisTotal = toc(tTransferRedis);
        fprintf(1, '\n  Total time to transfer to Redis  : %.3fs', tTransferRedisTotal);
        
        %% Delete redisRequestKey from redis
%         fprintf(1,'\n  Deleting key: %s',totalJobsKey)
%         [~, statusDel] = tools.redis.redisDEL(R,redisRequestKey);
%         if ~strcmp(statusDel,'OK')
%             ME = MException('eduTool:delKeyRedisFailure',...
%                 '%s : could not be deleted from redis',redisRequestKey);
%             throw(ME);
%         end
        
        %% Fill in the redisResponseKey with the redisResponseKeyValue
        fprintf(1,'\n  Setting key: %s',redisResponseKey)
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
        
        %% Flush data from input and output buffer of R
        flushinput(R)
        flushoutput(R)
        
        if breakLoop
            tools.redis.redisDisconnect(R);
            disp('HARD STOP OF CALCULATOR')
            break
        end

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
        
        % Flush data from input and output buffer of R
        flushinput(R)
        flushoutput(R)
        
        tools.redis.redisDisconnect(R);
        
        disp(errorMessage)
       
    end
end