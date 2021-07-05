function [updateStuct,statusSet] = updateJsonExperimentProgress(expControl,...
    type,messages,progress)

redisConnection = expControl.redis.R;

if nargin<=3
    type            = 'error';
    message         = '';
    messageAdmins   = '';
    progress        = '';
else
    message         = messages.message;
    messageAdmins   = messages.messageAdmins;
end

%% update the correct redis key 

updateStuct.type            = type;
updateStuct.message         = message;
updateStuct.messageAdmins   = messageAdmins;
updateStuct.progress        = progress;

if strcmp(type,'info')
    % if type = info, then update the 
    % EXPERIMENT_{InstanceID}_{ExperimentID}_INFO redis key 
    redisKey        = expControl.redis.keys.experimentInfoRedisKey;

    infoJson        = jsonencode(updateStuct);

    [~, statusSet]  = tools.redis.redisSetJson(redisConnection,redisKey,infoJson);
else
    % for all the other types, update the
    % EXPERIMENT_{InstanceID}_{ExperimentID}_UPDATES redis key 
    redisKey        = expControl.redis.keys.experimentUpdatesRedisKey;

    updateJson    	= jsonencode(updateStuct);
    
    [~, statusSet]  = tools.redis.redisSetJson(redisConnection,redisKey,updateJson);
end

%% update the EXPERIMENT_{InstanceID} redis key for certain statuses
% update the redis key when status = started, error, cancelled-error,
% confirm or finished
if ismember(type,{'started','error','cancelled-error','confirm','finished'})
    
    % Read the EXPERIMENT_{InstanceID} redis key, just in case there is a
    % change to its value. Get from redis the json file that holds all the 
    % experiments
    [experimentsStr,~,~]    = tools.redis.redisGet(...
        expControl.redis.R,...
        expControl.redis.keys.experimentsRedisKey);

    experimentArray = jsondecode(experimentsStr{1,1}); 
        
    latestExperimentID  = expControl.latestExperimentID;  
    idExp = find(strcmp({experimentArray.experiment_id},...
        num2str(latestExperimentID)));
    
    % Only if the experiment has not been cancelled by the user, update its
    % status
    if ~strcmp(experimentArray(idExp).status,'cancelled')    
        experimentArray(idExp).status = type; %#ok

        tools.redis.redisSetJsonWrapper(redisConnection,...
            expControl.redis.keys.experimentsRedisKey,experimentArray,1);
        if ismember(type,{'cancelled-error','confirm'})
            error('MESSAGE-BY-PLATFORM');
        end
    end
end