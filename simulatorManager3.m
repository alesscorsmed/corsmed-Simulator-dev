function simulatorManager3(redisAddress,redisPort,...
    redisRequestKey,redisResponseKey,redisUpdatesKey,uniqueID,jobID,...
    dummyData)
%
% SERVICES.SIMULATORMANAGER
%
%	Mock up of the SIMULATOR manager service
%
% INPUT
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'simulatorManager';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

if (nargin < 8)
    dummyData = 0;
end

useMatlabBuiltinRedisFunction = 0;

try

    %% Establish connection with Redis
    % Connect to a redis server
    redis.address   = redisAddress;
    redis.port      = redisPort;
    
    disp('Trying to connect to redis')
    redis.address       = redisAddress;
    redis.port          = redisPort;
    inputBufferSize     = 10000000;
    outputBufferSize    = 10000000;

    R = tools.redis.redisEstablishConnection(redis.address,redis.port,...
        inputBufferSize,outputBufferSize);
    
    %%
    % Check if the key exists and it is not empty
    fprintf(1,'\n  Getting data for key: %s',redisRequestKey)
    [redisRequest,~,redisReport] = tools.redis.redisGet(R,redisRequestKey);

    % If there is no key, or the key is empty, continue
    if isempty(redisRequest) || strcmp(redisReport,'ERROR - NONEXISTANT KEY')
        ME = MException('eduTool:emptySimReqRedisKey',...
            '%s : this key is either empty or does not exist',redisRequestKey);
        throw(ME);
    elseif iscell(redisRequest)
        if isempty(redisRequest{1,1})
            ME = MException('eduTool:emptySimReqRedisKey',...
            '%s : this key is either empty',redisRequestKey);
        throw(ME);
        end
    end
    
    tTotal = tic();
    %% GET DATA FROM REDIS
    fprintf(1, '\n Getting data from redis for %s experiment (worker %s)',...
        uniqueID,jobID);
    fprintf(1, '\n');
    
%     spinModelKey        = [uniqueID,'_spinModel'];
    pulseSequenceKey    = strcat(uniqueID,'_pulseSequence');
    motionModelKey      = strcat(uniqueID,'_motionModel');
    expControlKey       = strcat(uniqueID,'_expControl');
    simJobKey           = strcat(uniqueID,'_simJob',jobID);
    
    
    % Flush data from input and output buffer of R
    flushinput(R)
    flushoutput(R)

    % Get data from redis
    fprintf(1,'\n  Getting data for key: %s',pulseSequenceKey)
    [pulseSequence,R]   = tools.redis.redisGetJsonWrapper(R,pulseSequenceKey);
    fprintf(1,'\n  Getting data for key: %s',motionModelKey)
    [motionModel,R]     = tools.redis.redisGetJsonWrapper(R,motionModelKey);
    fprintf(1,'\n  Getting data for key: %s',expControlKey)
    [expControl,R]      = tools.redis.redisGetJsonWrapper(R,expControlKey);

    % Flush data from input and output buffer of R
    flushinput(R)
    flushoutput(R)

    if dummyData
        disp('SKIPPING REDIS: LOADING DATA FOR simJob FROM LOCAL SYSTEM')
        load('20210125_simJob')
    else
        fprintf(1,'\n  Getting data for key: %s',simJobKey)
        [simJob,R] = tools.redis.redisGetJsonWrapper(R,simJobKey);
    end
    
    % Construct spinModel TEMP @@@
    spinModel.numParts  = 1;
    spinModel.numJobs   = 1;
    spinModel.numSlices = 1; % @@@
    spinModel.data{1,1} = simJob;
    
    % Define the numJob that the current simulator instance will execute,
    % the R connection to redis and the redisUpdatesKey
    expControl.numJob           = str2num(jobID);
    expControl.R                = R;
    expControl.redisUpdatesKey  = redisUpdatesKey;
    
    % TEMP for addressing issues
    expControl.estimatedRunTime = seconds(0);

    %% RUN SIMULATOR
    
    fprintf(1, '\n Starting simulator service for %s experiment', uniqueID);
    fprintf(1, '\n');
    
    [simSignal] = eduTool.run.engineService( ...
        spinModel, pulseSequence, motionModel, expControl );

    %% SAVE DATA TO REDIS
    simSignalKey = strcat(uniqueID,'_simSignal_',jobID);    
    R = tools.redis.redisSetJsonWrapper(R,simSignalKey,simSignal);
            
    %% Delete redisRequestKey from redis
%     [~, statusDel] = tools.redis.redisDEL(R,redisRequestKey);
%     if ~strcmp(statusDel,'OK')
%         ME = MException('eduTool:delKeyRedisFailure',...
%             '%s : could not be deleted from redis',redisRequestKey);
%         throw(ME);
%     end

    %% REPORT OK
    fprintf(1, '\n');
    fprintf(1, '\n SIMULATOR MANAGER on %s DONE', redisRequestKey);
    fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
    fprintf(1, '\n');
    
    % Flush data from input and output buffer of R
    flushinput(R)
    flushoutput(R)
    
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
