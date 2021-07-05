function [acquisition] = initialize(expControl,sessionData)
%
% DATA.ACQUISITIONDATA.INITIALIZE
%
%	Function that initializes acquisition data structure.
%   Returns a acquisition with empty fields, or filed if connDB is passed.
%
%   This function is useful to define the fields.
%
% INPUT
%   expControl   struct with experiment info
%
% OUTPUT
%   acquisition   acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%
if nargin < 1
    fillData = 0;
else
    fillData = 1;
end

% generate empty basic data
acquisition.data            = data.acquisition.initializeData();
% main RF pulse
acquisition.mainRF          = data.acquisition.initializeRF(90);
% refocusing pulse, if required
acquisition.refRF           = data.acquisition.initializeRF(180);
% IR preparation info, if needed
acquisition.prepIR          = data.acquisition.initializeIR();
% PG encoding info, if needed
acquisition.encPG           = data.acquisition.initializePG();
% kSpaceInfo, information about the encodings
acquisition.kSpaceInfo  = [];
% reconstruction, information about the reconstruction (ISMRMRD format)
acquisition.reconstruction  = data.reconData.initialize();

% if expControl is passed, filled the ACQ data with the info
if fillData
    
    try 
        %% get main data from DB
        sqlQueryUser = ['SELECT pulseq_parameters.unique_name, ',...
            'pulseq_parameter_value.selected_value FROM `pulseq_parameter_value` ',...
            'RIGHT JOIN `pulseq_parameters` ON pulseq_parameter_value.parameter_id=',...
            'pulseq_parameters.id WHERE pulseq_parameter_value.pulseq_id=',...
            num2str(expControl.pulseqID)];
        sqlQueryResults = exec(expControl.connLocalDB, sqlQueryUser);
        bUser = fetch(sqlQueryResults);
        cellUser = bUser.Data; % this is a cell array. Convert it to a struct
        for i=1:size(cellUser,1)
            %disp(['structExper.',cellUser{i,1},'=',cellUser{i,2}]);
            structExper.(cellUser{i,1}) = cellUser{i,2};
        end
        % use remote data
        if isfield(sessionData, 'talkToRemoteDB') && sessionData.talkToRemoteDB
            sqlQueryAdminUserGroup = ['SELECT config FROM config_admin_usergroup',...
                ' WHERE usergroup=',num2str(sessionData.userGroupID)];
            sqlQueryResults = exec(sessionData.connRemoteDB, sqlQueryAdminUserGroup);
            bAdminUsergroup = fetch(sqlQueryResults);
            configAdminUsergroup  = bAdminUsergroup.Data{1,1};
            matches = strsplit(configAdminUsergroup,';');
            for i = 1:size(matches,2)-1
                splitOutput = strsplit(matches{1,i},'=');
                structExper.(splitOutput{1,1}) = splitOutput{1,2};
            end
            % LEVEL 3:
            sqlQueryAdminUser = ['SELECT config FROM config_admin_user',...
                ' WHERE user=',num2str(sessionData.userID)];
            sqlQueryResults = exec(sessionData.connRemoteDB, sqlQueryAdminUser);
            bAdminUser = fetch(sqlQueryResults);
            configAdminUser   = bAdminUser.Data{1,1};
            matches = strsplit(configAdminUser,';');
            for i = 1:size(matches,2)-1
                splitOutput = strsplit(matches{1,i},'=');
                structExper.(splitOutput{1,1}) = splitOutput{1,2};
            end
        end
        
        %% other info related to sequence
        sqlQuery = ['SELECT pulse_seq_scheme_id,seqnum FROM edt_tool_local.pulse_sequence',...
            ' WHERE exper_id=',num2str(expControl.experimentID)];
        sqlQueryResults = exec(expControl.connLocalDB, sqlQuery);
        sqlQueryResults = fetch(sqlQueryResults);
        pulseSeqSchem   = sqlQueryResults.Data{1,1};
        pulseSeqNum     = sqlQueryResults.Data{1,2};
        % 
        sqlQuery  = ['SELECT type, filename, analytical, is_3d',...
            ' FROM edt_tool_local.pulse_seq_scheme WHERE id =',...
            num2str(pulseSeqSchem)];
        sqlQueryResults     = exec(expControl.connLocalDB, sqlQuery);
        sqlQueryResults     = fetch(sqlQueryResults);
        pulseSeqType        = sqlQueryResults.Data{1,1};
        pulseSeqFamilyName  = sqlQueryResults.Data{1,2};
        pulseSeqAnalytical  = sqlQueryResults.Data{1,3};
        is3D                = sqlQueryResults.Data{1,4};
        
        %% query the points of the plane
        sqlQuery = ['SELECT sim_points,points ',...
            'FROM edt_tool_local.pulse_sequence WHERE exper_id=',...
            num2str(expControl.experimentID)];
        sqlQueryResults     = exec(expControl.connLocalDB, sqlQuery);
        sqlQueryResults     = fetch(sqlQueryResults);
        pointsAll           = sqlQueryResults.Data(:,1);
        pointsAll           = strsplit(pointsAll{1,1},'|');
        pointsAllFrontend   = sqlQueryResults.Data(:,2);
        pointsAllFrontend   = strsplit(pointsAllFrontend{1,1},'|');
        
    catch ME
        ME.identifier;
        ME.message;
        
        %% inform user
        msg = ['There was a connectivity error. ',...
            'The error has been reported for further review. ',...
            'Please try again and if error persist '...
            'contact your administrator.'];
        eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
            'cancelled-error',expControl.experimentID,expControl.pulseqID);
        
        %% send error in connection to DB for backend
        timeStamp = datestr(now,'yyyymmdd-HHMMSS');
        errorMessage = sprintf(['%s - User (%d) - Exper (%d) - ',...
            '(DB CONNECTION - acqData) Error in function %s() at line %d.',...
            '\n\nError Message:\n%s'], ...
            timeStamp,expControl.userID,expControl.experimentID,...
            ME.stack(1).name,ME.stack(1).line,ME.message);
        fprintf(1, '%s\n', errorMessage);
        
        % Write error message in db. Do not write this message if this
        % comes from CX, JV and GB instances
        if ~ismember(user_id,[790,933,1139])
            eduTool.frontend.notifyAdminForErrors(expControl.connLocalD,...
                errorMessage,expControl.instanceID)
        end
    end
    
    %% proceed to fill the data available
    
    % sequence data
    acquisition.data.pulseSeqSchem       = pulseSeqSchem;
    acquisition.data.pulseSeqNum         = pulseSeqNum;
    acquisition.data.pulseSeqFamilyName  = pulseSeqFamilyName;
    acquisition.data.pulseSeqType        = pulseSeqType;
    acquisition.data.pulseSeqAnalytical  = pulseSeqAnalytical;
    acquisition.data.is3D                = is3D;
    
    % plane
    acquisition.data.pointsAll           = pointsAll;
    acquisition.data.pointsAllFrontend   = pointsAllFrontend;
    
    % basic info, should always be available
    acquisition.data.numSlices      = str2num(structExper.slices);
    acquisition.data.sliceThickness = str2num(structExper.sliceThickness);
    if is3D && ~contains(lower(pulseSeqFamilyName), 'conc')
        acquisition.data.sliceGap = 0.0;
        acquisition.data.fovSE = acquisition.data.numSlices *...
            ( acquisition.data.sliceThickness + acquisition.data.sliceGap );
        acquisition.data.numSE = acquisition.data.numSlices;
        acquisition.data.matrixZ = acquisition.data.numSlices;
    else
        acquisition.data.sliceGap = str2num(structExper.slice_gap);
        acquisition.data.fovSE = acquisition.data.sliceThickness;
        acquisition.data.numSE = 0;
        acquisition.data.matrixZ = 1;
    end
    
    acquisition.data.matrixX = str2num(structExper.matrix_x);
    acquisition.data.matrixY = str2num(structExper.matrix_y);

    acquisition.data.fovFE = str2num(structExper.fov_x);
    acquisition.data.fovPE = str2num(structExper.fov_y);
    
    acquisition.data.numFE = acquisition.data.matrixX;
    acquisition.data.numPE = acquisition.data.matrixY;
    
    if isfield(structExper, 'te')
        acquisition.data.TE = str2num(structExper.te);
    end
    
    if isfield(structExper, 'tr')
        acquisition.data.TR = str2num(structExper.tr);
    end
    if isfield(structExper, 'mp_rage_tr')
        acquisition.data.shotTR = str2num(structExper.mp_rage_tr);
    else
        acquisition.data.shotTR = acquisition.data.TR;
    end
    
    acquisition.data.rxBW    = str2num(structExper.BW);
    if isfield(structExper, 'BWpx')
        acquisition.data.pixelBW = str2num(structExper.BWpx);
    end
    
    acquisition.data.NEX = str2num(structExper.NEX);
    
    acquisition.data.reconstructor = structExper.reconstructor;
    
    % main RF
    if isfield(structExper, 'RFflipAngle')
        acquisition.mainRF.flipAngle       = str2num(structExper.RFflipAngle);
    else
        acquisition.mainRF.flipAngle = 90;
    end
    acquisition.mainRF.duration        = str2num(structExper.RFduration);
    acquisition.mainRF.sliceThickness  = str2num(structExper.sliceThickness);
    % refocusing RF consistent if needed
    acquisition.refRF.duration        = str2num(structExper.RFduration);
    acquisition.refRF.sliceThickness  = str2num(structExper.sliceThickness);
    
    % optional entries
    if isfield(structExper, 'OOP_deltaTime')
        acquisition.data.numEchoes = 2;
        acquisition.data.deltaTimeOOP = str2num(structExper.OOP_deltaTime);
    end
    if isfield(structExper, 'kspaceshift')
        acquisition.data.kspaceshift  = structExper.kspaceshift;
    end
    if isfield(structExper, 'foldoversuppr')
        acquisition.data.foldoverSuppr  = structExper.foldoversuppr;
        if strcmpi(acquisition.data.foldoverSuppr, 'yes')
            acquisition.data.samplingFactorPE   = 2;
        else
            acquisition.data.samplingFactorPE   = 1;
        end
    end
    if isfield(structExper, 'foldoverdir')
        acquisition.data.foldoverDir    = structExper.foldoverdir;
    end
    if isfield(structExper, 'etl')
        acquisition.data.ETL = str2num(structExper.etl);
    end
    if isfield(structExper, 'concatenations')
        acquisition.data.concatenations = str2num(structExper.concatenations);
    end
    if isfield(structExper, 'parallelImaging')
        acquisition.data.parallelImaging = structExper.parallelImaging;
        if isfield(structExper, 'rfactor')
            acquisition.data.rFactor = str2num(structExper.rfactor);
        end
        if strcmpi(structExper.parallelImaging, 'sense')
            acquisition.data.reconstructor   = 'sense';
        end
    end
    if isfield(structExper, 'partialFourier')
        acquisition.data.partialFourier = structExper.partialFourier;
    end
    if isfield(structExper, 'partialFourierFactor')
        acquisition.data.fFactor = str2num(structExper.partialFourierFactor);
    end
    

    if isfield(structExper, 'refoc_pulse_angle')
        acquisition.refRF.flipAngle = str2num(structExper.refoc_pulse_angle);
    end
    if isfield(structExper, 'ir_time')
        acquisition.prepIR.Apply     = 1;
        acquisition.prepIR.TI        = str2num(structExper.ir_time);
    end
    if isfield(structExper, 'preparation_rf_angle')
        acquisition.prepIR.flipAngle = str2num(structExper.preparation_rf_angle);
    else
        acquisition.prepIR.flipAngle = 180; % in degrees
    end
    if isfield(structExper, 'ir_duration')    
        acquisition.prepIR.duration  = str2num(structExper.ir_duration);
    end
    if isfield(structExper, 'ir_type')
        acquisition.prepIR.type      = structExper.ir_type;
    end
    if isfield(structExper, 'fatsat')
        acquisition.data.fatsat  = structExper.fatsat;
    end

    % Fill in the reconstruction data
    acquisition.reconstruction.ismrmrd.header.encoding.encodedSpace.matrixSize.x  = ...
        acquisition.data.numFE;
    acquisition.reconstruction.ismrmrd.header.encoding.encodedSpace.matrixSize.y  = ...
        acquisition.data.numPE;
    acquisition.reconstruction.ismrmrd.header.encoding.encodedSpace.matrixSize.z  = ...
        acquisition.data.numSE;
    
    acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.x  = ...
        acquisition.data.matrixX;
    acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.y  = ...
        acquisition.data.matrixY;
    acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.z  = ...
        acquisition.data.matrixZ;

    acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.x = ...
        1000*acquisition.data.fovFE;
    acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.y = ...
        1000*acquisition.data.fovPE;
    acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.z = ...
        1000*acquisition.data.fovSE;

end

