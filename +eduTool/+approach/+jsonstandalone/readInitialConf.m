function [sessionData] = readInitialConf()
%
% EDUTOOL.APPROACH.JSONSTANDALONE.READINITIALCONF
%
%	reads tags and information from 
%
% INPUT
%
%   none
%
% OUTPUT
%   sessionData     session info and attributes from AWS instance
%
%========================  CORSMED AB © 2020 ==============================
%
%
functionName = 'eduTool.approach.jsonstandalone.readInitialConf';
% time it
tTotal = tic();
fprintf(1, '\n%s : start', functionName);
try
    %% initialize the session data structure
    [sessionData] = data.sessionData.initializeSession();
    %% get the info for instance id
    [~,instance_id] = system('sudo ec2metadata --instance-id');
    instance_id     = strtrim(instance_id);
    
    sessionData.instanceData = aws.ec2.jsonstandalone.readInstanceTags(instance_id);

    %% get the info for redis connection
    [sessionData] = eduTool.setup.connectRedis(sessionData);

    %% read the initial json file
    % redis key name for reading the session details
    redisRequestKey = ['SIMULATOR_',instance_id];
    
    % read redis key
    [instanceData,sessionData.redis.R,redisReport] = ...
        tools.redis.redisGetJsonWrapper(sessionData.redis.R,redisRequestKey);

    % If there is no key, or the key is empty, continue
    if isempty(instanceData) && ...
        (strcmp(redisReport,'ERROR - NONEXISTANT KEY') || ...
        strcmp(redisReport,'EMPTY KEY'))
        
        ME = MException('eduTool:readInitialRedisKeyFailure',...
            '%s : could not read this key from redis',redisRequestKey);
        throw(ME);
    end
    
    % redis key name for updating simulator status
    simulatorUpdatesRedisKey    = [redisRequestKey,'_UPDATES'];
    experimentsRedisKey         = ['EXPERIMENT_',instance_id];    
    
    sessionData.redis.keys.redisRequestKey          = redisRequestKey;
    sessionData.redis.keys.simulatorUpdatesRedisKey = simulatorUpdatesRedisKey;
    sessionData.redis.keys.experimentsRedisKey      = experimentsRedisKey;
    
    %% convert data and store in sessionData
    sessionData.instanceID     = instance_id;
    sessionData.userID         = str2num(instanceData.userID);
    
    sessionData.AWStagUserID   = 0;
    
    if isfield(instanceData,'Name')
        sessionData.instanceName    = instanceData.Name;
    else
        sessionData.instanceName    = 'N/A';
    end
    
    if isfield(instanceData,'courseID')
        sessionData.courseID        = str2num(instanceData.courseID);
    else
        sessionData.courseID        = 'N/A';
    end
    
    if isfield(instanceData,'anatomicalID')
        sessionData.anatomicalID    = str2num(instanceData.anatomicalID);
    else
        sessionData.anatomicalID    = 6;
    end

    %% drivers and dev info
    if isfield(instanceData,'CUDA')
        sessionData.cudaVersion     = str2num(instanceData.CUDA);
    else
        sessionData.cudaVersion     = 10;
    end
    if isfield(instanceData,'parfeval')
        sessionData.parfeval        = str2num(instanceData.parfeval);
    else
        sessionData.parfeval        = 1;
    end
    if isfield(instanceData,'PYTHON')
        sessionData.pythonVersion   = instanceData.PYTHON;
    else
        sessionData.pythonVersion   = '3.8';
    end
    if isfield(instanceData,'Development')
        sessionData.developmentUse  = str2num(instanceData.Development);
    else
        sessionData.developmentUse  = 0;
    end

    %% report
    fprintf(1, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(1, '\n  User ID    : %d', sessionData.userID);
    fprintf(1, '\n  Instance   : %s (%s)', sessionData.instanceID, sessionData.instanceName);
    fprintf(1, '\n  Course ID  : %d', sessionData.courseID);
    fprintf(1, '\n');
    fprintf(1, '\n');

catch ME
    
    ME.identifier;
    ME.message;
    
    %% send error in connection to DB for backend
    errorMessage = tools.printErrorMessage(sessionData,ME);
end