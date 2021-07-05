function [sessionData, tagsStruct] = initializeInstance(appSpecs)
%
% EDUTOOL.INITIALIZEINSTANCE
%
%	Initializes the instace, connecting to DB 
%   and retrieving session data
%
% INPUT
%
%   none
%
% OUTPUT
%
%   sessionData    struct with session and instance info
%   tagsStruct     struct with necessary data from AWS tags
%
%========================  CORSMED AB © 2020 ==============================
%

functionName = 'eduTool.initializeInstance';
% time it
tTotal = tic();
fprintf(1, '\n%s : start', functionName);

APP = appSpecs.application;
VERSION = appSpecs.version;
if ~isfield(appSpecs,'approach')
    APPROACH = '';
else
    APPROACH = appSpecs.approach;
end
if ~isfield(appSpecs,'mode')
    MODE = '';
else
    MODE = appSpecs.mode;
end

% try
if strcmp(APPROACH,'jsonstandalone')
    sessionData = eduTool.approach.jsonstandalone.readInitialConf();
    
    % assign application and version
    sessionData.application	= APP;
    sessionData.versionNum  = VERSION;
    sessionData.approach    = APPROACH;
    sessionData.mode        = MODE;  
    
    instanceID = sessionData.instanceID;
    
    tagsStruct = '';
    
    % Update the folderSystem
    jsonPath        = ['/efs-mount-point/S20/INPUTS/configuration/',...
        'edutool-jsonstandalone-folderSystem.json'];
    fid             = fopen(jsonPath);
    rawJson         = fread(fid,inf);
    strJson         = char(rawJson');
    fclose(fid);

    sessionData.folderSystem  = jsondecode(strJson);
elseif strcmp(APPROACH,'centralized')
    [~,instance_id] = system('sudo ec2metadata --instance-id');
    instance_id     = strtrim(instance_id);
    
    %% prepare initial setup
    sessionData.instanceID      = instance_id;
    sessionData.instanceName    = '';
    sessionData.courseID        = appSpecs.courseID;
    sessionData.userID          = appSpecs.userID;
    sessionData.anatomicalID    = appSpecs.anatomicalID;
    sessionData.AWStagUserID    = 0;

    %% drivers and dev info
    sessionData.cudaVersion     = 10;
    sessionData.parfeval        = 1;
    sessionData.pythonVersion   = '3.8';
    sessionData.developmentUse  = 1;
    
    %% dummy data for tagStruct
    tagsStruct = '';
    
    %% connect to DB
    [sessionData] = eduTool.setup.connectDatabase(sessionData);
    
    %% Update DB that backend has started
    [instanceID] = eduTool.frontend.updateBackendStartCounts(sessionData.connLocalDB);
    
    %% assign application and version
    sessionData.application	= APP;
    sessionData.versionNum  = VERSION;
    sessionData.approach    = APPROACH;
    sessionData.mode        = MODE;    
    
    %% NEW FOLDER SYSTEM: structured
    sessionData.folderSystem.baseFolder = '/efs-mount-point/S20';
    % from the base folder, generate the different INPUT folders
    sessionData.folderSystem.anatomicalModelFolder = ...
        sprintf('%s/INPUTS/anatomical/%s/', ...
        sessionData.folderSystem.baseFolder, sessionData.application);
    sessionData.folderSystem.coilModelFolder = ...
        sprintf('%s/INPUTS/coils/%s/', ...
        sessionData.folderSystem.baseFolder, sessionData.application);
    sessionData.folderSystem.mrSystemModelFolder = ...
        sprintf('%s/INPUTS/system/%s/', ...
        sessionData.folderSystem.baseFolder, sessionData.application);
else
    %% initialize AWS instance & Session attributes
    [sessionData,tagsStruct] = aws.ec2.readInstanceTags();
    %% connect to DB
    [sessionData] = eduTool.setup.connectDatabase(sessionData,tagsStruct);
    
    %% Update DB that backend has started
    [instanceID] = eduTool.frontend.updateBackendStartCounts(sessionData.connLocalDB);
    
    %% assign application and version
    sessionData.application	= APP;
    sessionData.versionNum  = VERSION;
    sessionData.approach    = APPROACH;
    sessionData.mode        = MODE;    
    
    %% NEW FOLDER SYSTEM: structured
    sessionData.folderSystem.baseFolder = '/efs-mount-point/S20';
    % from the base folder, generate the different INPUT folders
    sessionData.folderSystem.anatomicalModelFolder = ...
        sprintf('%s/INPUTS/anatomical/%s/', ...
        sessionData.folderSystem.baseFolder, sessionData.application);
    sessionData.folderSystem.coilModelFolder = ...
        sprintf('%s/INPUTS/coils/%s/', ...
        sessionData.folderSystem.baseFolder, sessionData.application);
    sessionData.folderSystem.mrSystemModelFolder = ...
        sprintf('%s/INPUTS/system/%s/', ...
        sessionData.folderSystem.baseFolder, sessionData.application);
end
     
% catch
%     ME = MException('EduTool:instanceError',...
%         '%s : unable to initialize the instance',functionName);
%     throw(ME);
% end

%% report
fprintf(1, ...
    '\n%s : done, elapsed time %.3fs',...
    functionName, toc(tTotal));
fprintf(1, '\n  INSTANCE %s READY', instanceID);
fprintf(1, '\n');
