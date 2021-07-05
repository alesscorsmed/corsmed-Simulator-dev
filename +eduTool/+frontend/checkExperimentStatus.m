function checkExperimentStatus(expControl,experimentID)

%FRONTEND.CHECKEXPERIMENTSTATUS
%
% checks the possibility of cancelled experiment
%
% INPUT
%   experimentID        ID of the corresponding experiment
%   expControl          the structure that holds the DB connection object
%
% OUTPUT
%   none
%
% ========================= CORSMED AB @ 2020 ============================
%

if isfield(expControl, 'connLocalDB') ...
        && ~isempty(expControl.connLocalDB)
    sqlquery = ['SELECT status FROM',...
        ' edt_tool_local.experiments WHERE id=',num2str(experimentID)];
    sqlquery_results    = exec(expControl.connLocalDB, sqlquery);
    sqlquery_results    = fetch(sqlquery_results);
    status              = sqlquery_results.Data{1,1};
    
    if strcmp(status,'cancelled')
        statusText = 'The virtual MR scanner is now ready for the next experiment';
        eduTool.frontend.updateScannerStatus(expControl,statusText);
        msg = 'Experiment cancelled by the user.';
        error('CANCELLED-BY-USER');
    end
    
elseif isfield(expControl,'approach') && strcmp(expControl.approach,'jsonstandalone')
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
    if strcmp(experimentArray(idExp).status,'cancelled')
        error('CANCELLED-BY-USER');
    end    
end