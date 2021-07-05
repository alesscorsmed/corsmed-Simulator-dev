function [sessionData, tagsStruct] = initializeLocalInstance(APP, VERSION)
%
% EDUTOOL.INITIALIZELOCALINSTANCE
%
%	Initializes the instace, with no DB connection 
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
functionName = 'eduTool.initializeLocalInstance';
% time it
tTotal = tic();
fprintf(1, '\n%s : start', functionName);

% try
    %% initialize AWS instance & Session attributes
    [sessionData,tagsStruct] = aws.ec2.readInstanceTags();
    
    %% assign application and version
    sessionData.application	= APP;
    sessionData.versionNum  = VERSION;
    
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
     
% catch
%     ME = MException('EduTool:instanceError',...
%         '%s : unable to initialize the instance',functionName);
%     throw(ME);
% end

%% report
fprintf(1, ...
    '\n%s : done, elapsed time %.3fs',...
    functionName, toc(tTotal));
fprintf(1, '\n  INSTANCE %s READY', sessionData.instanceID);
fprintf(1, '\n');
