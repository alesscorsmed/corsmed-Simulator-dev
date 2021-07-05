function [sessionData,tagsStruct] = readInstanceTags()
%
% EDUTOOL.AWS.EC2.READINSTANCETAGS
%
%	reads tags and information from the AWS instance.
%
% INPUT
%
%   none
%
% OUTPUT
%   tagsStruct      struct with necessary data from AWS tags
%   sessionData     session info and attributes from AWS instance
%
%========================  CORSMED AB © 2020 ==============================
%
%
functionName = 'eduTool.aws.ec2.readInstanceTags';
% time it
tTotal = tic();
fprintf(1, '\n%s : start', functionName);
%% initialize the session data structure
[sessionData] = data.sessionData.initializeSession();
%% get the info
[~,instance_id] = system('sudo ec2metadata --instance-id');
instance_id     = strtrim(instance_id);
% Read the assigned tags
command2        = ['sudo aws ec2 describe-instances --instance-ids ',...
    instance_id,' --query "Reservations[].Instances[].Tags"'];
[~,tags]        = system(command2);
% decode
tagsStruct      = jsondecode(tags);
%% assign info to structure
for i=1:size(tagsStruct,2)
    strName = regexprep(tagsStruct(i).Key,'-','');
    instanceData.(strName) = tagsStruct(i).Value;
end

%% convert data and store in sessionData
sessionData.instanceID     = instance_id;
sessionData.instanceName   = instanceData.Name;
sessionData.courseID       = str2num(instanceData.Courseid);
sessionData.userID         = str2num(instanceData.Userid);
sessionData.AWStagUserID   = 0;

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