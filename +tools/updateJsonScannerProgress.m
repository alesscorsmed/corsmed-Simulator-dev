function [updateStuct,statusSet] = updateJsonScannerProgress(expControl,...
    type,message,progress)

redisConnection = expControl.redis.R;

if nargin==4
    progress = '';
elseif nargin<=3
    type        = 'error';
    message     = '';
    progress    = '';
end

%% update the correct redis key

% for all the other types, update the
% SIMULATOR_{InstanceID}_UPDATES redis key 
redisKey        = expControl.redis.keys.simulatorUpdatesRedisKey;

updateStuct.type        = type;
updateStuct.message     = message;
updateStuct.progress    = progress;

updateJson  = jsonencode(updateStuct);

[~, statusSet] = tools.redis.redisSetJson(redisConnection,redisKey,updateJson);