function [anatomicalModel,coilSystem,sessionData] = loadCourseData(...
    sessionData,tagsStruct)
%
% BACKEND.LOADCOURSEDATA
%
%	Function that initializes the data for a course.
%
% INPUT
%   courseID   course id
%
% OUTPUT
%   anatomicalModel   structure with model data
%   coilSystem        structure with coil models
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.setup.loadCourseData';
if (nargin < 1)
    anatomicalID                = 6;
    sessionData.anatomicalID    = anatomicalID;
elseif nargin == 1 || isempty(tagsStruct)
    anatomicalID                = sessionData.anatomicalID;
else
    struct_index                = find(strcmp({tagsStruct.Key}, 'Anatomical-id')==1);
    anatomicalID                = str2double(tagsStruct(struct_index).Value);
    sessionData.anatomicalID    = anatomicalID;
end

%% report start
fid = 1;
tTotal = tic();
fprintf(fid, '\n\n%s : start', functionName);

try
    
    %% load anatomical model
    statusMessage = 'Loading the anatomical model...';
    eduTool.frontend.updateScannerStatus(sessionData,statusMessage);
        
    [anatomicalModel] = data.models.initializeAnatomical(anatomicalID,...
        sessionData.folderSystem.anatomicalModelFolder, ...
        sessionData.application);
    
    %% load coil system
    eduTool.frontend.updateScannerStatus(sessionData, ...
        'Coil calibration and setup...');
        
    [coilSystem] = data.models.initializeCoils(anatomicalID,...
        sessionData.folderSystem.coilModelFolder);
    
    %% precompute SAR for tx coils
    eduTool.frontend.updateScannerStatus(sessionData, ...
        'MR safety checks...');

    [coilSystem] = coils.precomputeMRSafety(coilSystem, anatomicalModel);
    %% report
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for course ID %d',...
        functionName, anatomicalID );
    fprintf(fid, '\n  Scanner ready' );
    fprintf(fid, '\n    Active coil     %s', coilSystem.activeRx);
    fprintf(fid, '\n    Patient name    %s', anatomicalModel.name);
    fprintf(fid, '\n  Loading Time      %.3fs', tTotal);
    fprintf(fid, '\n\n');
    
catch ME
    
    ME.identifier;
    ME.message;
    
    %% send error in connection to DB for backend
    errorMessage = tools.printErrorMessage(sessionData,ME);
    
    % Write error message in db. Do not write this message if this
    % comes from CX, JV and GB instances
    if ~ismember(sessionData.userID,[790,933,1139]) && updateStatusThroughRedis==0
        eduTool.frontend.notifyAdminForErrors(sessionData.connLocalDB,...
            errorMessage,sessionData.instanceID)
    end
end