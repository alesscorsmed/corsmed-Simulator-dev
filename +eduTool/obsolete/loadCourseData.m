function [anatomicalModel,mrSystem,coilSystem] = loadCourseData(...
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
    anatomicalID = 6;
elseif nargin == 1
    anatomicalID = sessionData.anatomicalID;
else
    struct_index    = find(strcmp({tagsStruct.Key}, 'Anatomical-id')==1);
    anatomicalID 	= str2double(tagsStruct(struct_index).Value);
end

%% report start
fid = 1;
tTotal = tic();
fprintf(fid, '\n\n%s : start', functionName);

try
    %% initialize MR system
    if isfield(sessionData,'connLocalDB')
        statusMessage = 'Booting up the MR scanner...';
        eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
            statusMessage);
    end
    [mrSystem] = data.mrSystem.initialize();
    
    %% load anatomical model
    if isfield(sessionData,'connLocalDB')
        statusMessage = 'Loading the anatomical model...';
        eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
            statusMessage);
    end
    [anatomicalModel] = data.anatomicalModel.initialize(anatomicalID,...
        sessionData.folderSystem.anatomicalModelFolder, ...
        sessionData.application);
    
    %% load coil system
    if isfield(sessionData,'connLocalDB')
        eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
            'Coil calibration and setup...');
    end
    [coilSystem] = data.coilModel.initialize(anatomicalID,...
        sessionData.folderSystem.coilModelFolder);
    
    %% precompute SAR for tx coils
    if isfield(sessionData,'connLocalDB')
        eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
            'MR safety checks...');
    end
    [coilSystem] = coils.precomputeMRSafety(coilSystem, anatomicalModel);
    
    %% report
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for course ID %d',...
        functionName, anatomicalID );
    fprintf(fid, '\n  Scanner ready' );
    fprintf(fid, '\n    Scanner model   %s', mrSystem.name );
    fprintf(fid, '\n    Active coil     %s', coilSystem.activeRx);
    fprintf(fid, '\n    Patient name    %s', anatomicalModel.name);
    fprintf(fid, '\n  Booting Time      %.3fs', tTotal);
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