function [expControl] = dbQueryExpControl(expControl,sessionData)
%
% EDUTOOL.SETUP.DBQUERYEXPCONTROL
%
%	Function that loads expControl struct from DB queries
%
% INPUT
%   expControl   structure with acquisition parameters
%   sessionData  struct with the sessionData
%
% OUTPUT
%   expControl   structure with acquisition parameters
%
%========================  CORSMED AB © 2020 ==============================
%

try
    %% get the experiment data
    experiment = eduTool.frontend.expectExperiment(sessionData.connLocalDB);
    % if empty, no experiment yet, return empty expControl
    if strcmp(experiment.Data,'No Data')
        expControl = [];
        return;
    end
    
    %% session data
    expControl.application      = 'edutool';
    expControl.versionNum       = sessionData.versionNum;
    expControl.instanceID       = sessionData.instanceID;
    expControl.courseID         = sessionData.courseID;
    expControl.connLocalDB      = sessionData.connLocalDB;
    expControl.anatomicalID     = sessionData.anatomicalID;
    
    %% get experiment data
    expControl.experimentID    = experiment.Data{1,1};
    expControl.remotedbID      = experiment.Data{1,2};
    expControl.status          = experiment.Data{1,3};
    expControl.experInfo       = experiment.Data{1,4};
    expControl.reconstructor   = experiment.Data{1,5};
    expControl.reconInfo       = experiment.Data{1,6};
    expControl.pulseqMatFile   = experiment.Data{1,7};
    expControl.pulseqID        = experiment.Data{1,8};
    
    if ~isempty(expControl.pulseqMatFile)
        expControl.pathPulseSeq = ...
            [sessionData.folderSystem.pulseSequenceFolder,...
            filesep, expControl.pulseqMatFile];
    else
        expControl.pathPulseSeq = '';
    end
    
    %% user ID
    sqlquery = ['SELECT user_id FROM',...
        ' edt_tool_local.experiments WHERE id=',num2str(expControl.experimentID)];
    sqlqueryResults = exec(expControl.connLocalDB, sqlquery);
    b = fetch(sqlqueryResults);
    expControl.userID  = b.Data{1,1};
    
    %% debugging and development mode
    % comes from CX, JV and GB instances
    if ismember(sessionData.userID,[790,933,1139])
        expControl.debug.devMode      = 1;
        expControl.debug.debugMode    = 1;
    else
        expControl.debug.devMode      = 0;
        expControl.debug.debugMode    = 0;
    end
    expControl.debug.debugFile    = [];
    % expControl.debug.formatTime   = 'yyyymmdd-HHMMSS';
    % expControl.debug.timeStamp    = datestr(now,'yyyymmdd-HHMMSS');
    
    %% folder system
    expControl.folderSystem = sessionData.folderSystem;
    % update kernel path
    expControl.simulation.kernelFolder = ...
        sprintf('%s/PROJECTS/%s/kernels/',...
        expControl.folderSystem.baseFolder, expControl.application);
    % default kernel
    expControl.simulation.kernelPtx = sprintf('%s%s.ptx', ...
        expControl.simulation.kernelFolder, expControl.simulation.kernelVer);
    % upgrade with more info
    expControl.folderSystem.errorFolder = ...
        sprintf('%s/ERRORS/%s/', ...
        expControl.folderSystem.baseFolder, expControl.application);
    expControl.folderSystem.statsFolder = ...
        sprintf('%s/STATS/%s/', ...
        expControl.folderSystem.baseFolder, expControl.application);
    expControl.folderSystem.experimentFolder = ...
        sprintf('%s/RESULTS/%s/user_%d/', ...
        expControl.folderSystem.baseFolder, ...
        expControl.application, expControl.userID);
    if ~isfolder(expControl.folderSystem.experimentFolder)
        system(sprintf('mkdir %s',expControl.folderSystem.experimentFolder));
    end
    
    % Creation of new empty Struct
    for i=1:4
        for j=1:4
            expControl.stack_struct.slice.image.corners(i).line(j).text = "";
        end
    end
    
    sqlquery2           = ['SELECT selected_value FROM',...
        ' edt_tool_local.global_configuration WHERE name=''outer_fov_ratio'''];
    sqlqueryResults2   = exec(expControl.connLocalDB, sqlquery2);
    sqlqueryResults2   = fetch(sqlqueryResults2);
    
    expControl.outerFOVratio = str2num(sqlqueryResults2.Data{1,1});
    
    
    %% get simulation data
    sqlQuery        = ['SELECT selected_value FROM ',...
        'edt_tool_local.global_configuration WHERE',...
        ' name IN (''performance_gridZ_sliceThickness'',',...
        '''PDinhomogeneity'',''isotropic_grid'',''analytical_simulation'',',...
        '''fast_algorithm'',''coil'',''isocenter'',''motion'',''motion_rot_angle'',',...
        '''motion_rot_freq'',''motion_trans_magn'',''motion_trans_freq'',',...
        '''motion_trans_axis'',''kernel_threads'',''kernel_blocks'',',...
        '''gridSizeOption'',''activate_cs'',''dwell_time'',',...
        '''coil_basic_or_advanced'',''deactivateGx'',''deactivateGy'',',...
        '''deactivateGz'',''simulationMode'',''simulationKernel'',',...
        '''advanced_notifications'',''pulse_seq_generator'',''noise'') ORDER BY FIELD ',...
        '(name,''performance_gridZ_sliceThickness'',''PDinhomogeneity'',',...
        '''isotropic_grid'',''analytical_simulation'',''fast_algorithm'',',...
        '''coil'',''isocenter'',''motion'',''motion_rot_angle'',',...
        '''motion_rot_freq'',''motion_trans_magn'',''motion_trans_freq'',',...
        '''motion_trans_axis'',''kernel_threads'',''kernel_blocks'',',...
        '''gridSizeOption'',''activate_cs'',''dwell_time'',',...
        '''coil_basic_or_advanced'',''deactivateGx'',''deactivateGy'',',...
        '''deactivateGz'',''simulationMode'',''simulationKernel'',',...
        '''advanced_notifications'',''pulse_seq_generator'',''noise'')'];
    % query
    sqlQueryResults    = exec(sessionData.connLocalDB, sqlQuery);
    perfTypeInfo       = fetch(sqlQueryResults);
    
    %% update expControl (notice numbers of Data not in order!)
        
    % simulation related
    expControl.simulation.analyticalSim     = str2double(perfTypeInfo.Data{4,1}); % full analytical
    expControl.simulation.simulationEngine  = perfTypeInfo.Data{23,1}; % numerical / analytical
    expControl.simulation.simulationKernel  = perfTypeInfo.Data{24,1}; % latest / stable
    
    expControl.simulation.threads           = str2double(perfTypeInfo.Data{14,1});
    expControl.simulation.blocks            = str2double(perfTypeInfo.Data{15,1});
    expControl.simulation.activateCS        = str2double(perfTypeInfo.Data{17,1});
    
    % model related
    expControl.model.useSliceThickness  = str2double(perfTypeInfo.Data{1,1});
    expControl.model.pdInhomogeneity    = str2double(perfTypeInfo.Data{2,1});
    expControl.model.isotropicGrid      = str2double(perfTypeInfo.Data{3,1});
    expControl.model.coilType           = perfTypeInfo.Data{6,1};
    expControl.model.zIsocenter         = str2double(perfTypeInfo.Data{7,1})/100;
    if (expControl.anatomicalID == 6)
        % MIDA, fix isocenter, need to correct in Front End
        expControl.model.zIsocenter = 0.12;
    end
    expControl.model.coilMode           = perfTypeInfo.Data{19,1};  % It refers to the option coil_basic_or_advanced
    expControl.model.gridSizeOption     = perfTypeInfo.Data{16,1};
    switch lower(expControl.model.gridSizeOption)
        case '1mm'
            expControl.model.gridStep   = [1.0, 1.0, 1.0]*1e-3;
        case '0p5mm'
            expControl.model.gridStep   = [0.5, 0.5, 0.5]*1e-3;
        otherwise
            expControl.model.gridSizeOption = 'optimal';
            expControl.model.gridStep   = [0.5, 0.5, 1.0]*1e-3;
            expControl.simulation.simulationEngine = 'phasor';
    end
    % for native model resolution, use phasor
    if ~expControl.model.isotropicGrid
        expControl.simulation.simulationEngine = 'phasor';
    end
    
    % noise
    expControl.simulation.applyNoise        = str2num(perfTypeInfo.Data{27,1});
    
    % sequence related
    expControl.sequence.timeCompression = str2double(perfTypeInfo.Data{5,1}); % former fast algorithm
    expControl.sequence.dwellTime       = perfTypeInfo.Data{18,1};
    expControl.sequence.deactivateGx    = str2double(perfTypeInfo.Data{20,1});
    expControl.sequence.deactivateGy    = str2double(perfTypeInfo.Data{21,1});
    expControl.sequence.deactivateGz    = str2double(perfTypeInfo.Data{22,1});
    
    % motion related
    expControl.motionSpecs.pattern      = perfTypeInfo.Data{8,1};
    expControl.motionSpecs.rotAngle     = str2double(perfTypeInfo.Data{9,1});
    expControl.motionSpecs.rotFreq      = str2double(perfTypeInfo.Data{10,1});
    expControl.motionSpecs.transMag     = str2double(perfTypeInfo.Data{11,1});
    expControl.motionSpecs.transFreq    = str2double(perfTypeInfo.Data{12,1});
    expControl.motionSpecs.transAxis    = str2double(perfTypeInfo.Data{13,1});
    
    % standard and other
    expControl.advancedNotifications   = str2double(perfTypeInfo.Data{25,1});
    
    % Use the old (v1) or new (v2) pulse sequence generator
    pulseSeqGenerator = perfTypeInfo.Data{26,1};
    if strcmp(pulseSeqGenerator,'v1')
        expControl.useOldSequence = 1;
    else
        expControl.useOldSequence = 0;
    end

    
    
    % If anatomicalID = 9, get from db the default cardiac phase
    if (expControl.anatomicalID == 9)
        sqlQueryPhase = ['SELECT selected_value FROM ',...
            'edt_tool_local.global_configuration WHERE',...
            ' name IN (''defaultCardiacPhase'') ORDER BY FIELD ',...
            '(name,''defaultCardiacPhase'')'];
        sqlQueryPhaseResults = exec(sessionData.connLocalDB, sqlQueryPhase);
        perfTypeInfoPhase    = fetch(sqlQueryPhaseResults);

        expControl.model.cardiacPhase  = str2double(perfTypeInfoPhase.Data{1,1});
    end
    
    % Get the value for the menu item "Activate susceptibility"
    sqlQuerySusc = ['SELECT selected_value FROM ',...
        'edt_tool_local.global_configuration WHERE',...
        ' name IN (''susceptibility'') ORDER BY FIELD ',...
        '(name,''susceptibility'')'];
    sqlQuerySuscResults = exec(sessionData.connLocalDB, sqlQuerySusc);
    perfTypeInfoSusc    = fetch(sqlQuerySuscResults);
    suscData            = perfTypeInfoSusc.Data;
    if ~strcmp(suscData,'No Data')
        expControl.simulation.activateSusc = str2double(perfTypeInfoSusc.Data{1,1});
    end
    
catch ME
    ME.identifier;
    ME.message;
    
    %% inform user
    msg = ['There was a connectivity error. ',...
        'The error has been reported for further review. ',...
        'Please try again and if error persist '...
        'contact your administrator.'];
    eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
    
    %% send error in connection to DB for backend
    errorMessage = sprintf(['%s - User (%d) - Exper (N/A) - ',...
        '(DB CONNECTION - expControl) Error in function %s() at line %d.',...
        '\n\nError Message:\n%s'], ...
        expControl.timeStamp,expControl.userID,...
        ME.stack(1).name,ME.stack(1).line,ME.message);
    fprintf(1, '%s\n', errorMessage);
    
    % Write error message in db. Do not write this message if this
    % comes from CX, JV and GB instances
    if ~ismember(expControl.userID,[790,933,1139])
        eduTool.frontend.notifyAdminForErrors(sessionData.connLocalDB,...
            errorMessage,sessionData.instanceID)
    end
    
end

