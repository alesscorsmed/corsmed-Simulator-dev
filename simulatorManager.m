function simulatorManager(redisAddress,redisPort,...
    redisRequestKey,redisResponseKey,redisUpdatesKey,uniqueID,jobID)
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

try

    %% Establish connection with Redis
    % Connect to a redis server
    redis.address   = redisAddress;
    redis.port      = redisPort;

    R = tools.redis.redisEstablishConnection(redis.address,redis.port);
    
    %%
    % Check if the key exists and it is not empty
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
    
    % Establish redis connection using the MATLAB built-in function
    c = tools.redisMatlab.redisEstablishConnectionMatlab(...
        redis.address,redis.port);
        
    % Add spinModel in redis and verify that it is available in redis
    spinModelKey        = [uniqueID,'_spinModel'];
    pulseSequenceKey    = [uniqueID,'_pulseSequence'];
    motionModelKey      = [uniqueID,'_motionModel'];
%     acquisitionKey      = [uniqueID,'_acquisition'];
    expControlKey       = [uniqueID,'_expControl'];
    simJobKey           = [uniqueID,'_simJob',jobID];

    v = get(c,{spinModelKey,pulseSequenceKey,motionModelKey,...
        expControlKey,simJobKey});

    spinModel       = v{1,1};
    pulseSequence   = v{1,2};
    motionModel     = v{1,3};
    expControl      = v{1,4};
    simJob          = v{1,5};
    
    % Define the numJob that the current simulator instance will execute,
    % the R connection to redis and the redisUpdatesKey
    expControl.numJob           = str2num(jobID);
    expControl.R                = R;
    expControl.redisUpdatesKey  = redisUpdatesKey;

    %%
    
    fprintf(1, '\n Starting simulator service for %s experiment', uniqueID);
    fprintf(1, '\n');
    
    %%
    [simSignal] = eduTool.run.engineService( ...
        spinModel, pulseSequence, motionModel, expControl );

    %% SAVE DATA TO REDIS
    % Establish redis connection using the MATLAB built-in function
    c = tools.redisMatlab.redisEstablishConnectionMatlab(...
        redis.address,redis.port);

    % Add spinModel in redis and verify that it is available in redis
    simSignalKey = [uniqueID,'_simSignal_',jobID];
    put(c,simSignalKey,simSignal);
    if ~isKey(c,simSignalKey)
        ME = MException('eduTool:spinModelRedisSetFailure',...
            '%s : could not be stored in redis',simSignalKey);
        throw(ME);
    end
    %% Delete redisRequestKey from redis
    [~, statusDel] = tools.redis.redisDEL(R,redisRequestKey);
    if ~strcmp(statusDel,'OK')
        ME = MException('eduTool:delKeyRedisFailure',...
            '%s : could not be deleted from redis',redisRequestKey);
        throw(ME);
    end

    %% REPORT OK
    fprintf(1, '\n');
    fprintf(1, '\n SIMULATOR MANAGER on %s DONE', redisRequestKey);
    fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
    fprintf(1, '\n');
    
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
