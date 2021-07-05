function [expControl] = initialize(experiment,sessionData,application)
%
% DATA.EXPCONTROL.INITIALIZE
%
%	Function that initializes the expControl data structure,
%   with information for the experiment control.
%   Returns a simulationControl with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   sessionData   struct with sessionData
%
% OUTPUT
%   expControl   expControl structure
%
%========================  CORSMED AB © 2020 ==============================
%

if nargin < 1
    fillData = 0;
else
    fillData = 1;
end

timeStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
%% default data
expControl.name         = 'experiment';
expControl.userID       = '';
expControl.instanceID   = '';
expControl.experimentID = '';
expControl.courseID     = '';
expControl.pulseqID     = '';
expControl.versionNum   = '';
expControl.connLocalDB  = [];
expControl.progress     = 0;
expControl.folderSystem = [];
expControl.timeStamp    = timeStamp;
% change time stamp into file format yyyymmddHHMMSS
timeStamp = strrep(timeStamp,'-','');
timeStamp = strrep(timeStamp,':','');
timeStamp = strrep(timeStamp,' ','');
expControl.fileTimeStamp = timeStamp;

% debugging and development mode
expControl.debug.devMode      = 1;
expControl.debug.debugMode    = 1;
expControl.debug.debugFile    = [];
expControl.debug.waitBarBE    = 0;

% model related
expControl.model.useSliceThickness  = 1;
expControl.model.pdInhomogeneity    = 1;
expControl.model.isotropicGrid      = 1;
expControl.model.coilType           = 'optimal';
expControl.model.gridSizeOption     = '1mm';
expControl.model.gridStep           = [1.0, 1.0, 1.0]*1e-3;
expControl.model.coilMode           = 'basic';  % It refers to the option coil_basic_or_advanced
expControl.model.xExtendPct         = 0.1; % extend percentage in x in model interpolation (10% default)
expControl.model.rxCoilName         = 'none'; % no coil selected yet
expControl.model.zIsocenter         = 0.0; % isocenter for the coil positioning

% simulation related
expControl.simulation.numberOfSim       = 0; % number of simulations
expControl.simulation.analyticalSim     = 0; % full analytical
expControl.simulation.simulationEngine  = 'analytical'; % analytical / phasor / diffusion / numerical
expControl.simulation.precision         = 'single'; % single / double
expControl.simulation.odeMethod         = 'adaptiveExp'; % explicit / adaptiveExp / implicit / adaptiveImp
expControl.simulation.simulationKernel  = 'latest'; % latest / stable
expControl.simulation.activateCS        = 0;

% parallelization and GPU mode
expControl.simulation.gpuName       = 'None';
expControl.simulation.numGPU    	= 1;
expControl.simulation.threads     	= 256;
expControl.simulation.blocks     	= 256;
expControl.simulation.kernelVer  	= 'v26';
expControl.simulation.architecture  = 'sm70';
expControl.simulation.kernelFolder  = '/efs-mount-point/S20/PROJECTS/edutool/kernels/';
expControl.simulation.kernelPtx  	= sprintf('%s%s_%s.ptx', ...
    expControl.simulation.kernelFolder, expControl.simulation.kernelVer,...
    expControl.simulation.architecture);

% sequence related
expControl.sequence.dtGR            = 10e-6; % default time discretization for GR
expControl.sequence.dtRF            = 1e-6; % default time discretization for RF
expControl.sequence.tGuardRF        = 10e-6; % time guard to avoid RF/GR overlap
expControl.sequence.dwellTime       = 'dynamic';
expControl.sequence.timeCompression = 0; % former fast algorithm
expControl.sequence.minContextExc   = 0; % minimizes context change in kernel exec
expControl.sequence.deactivateGx    = 0;
expControl.sequence.deactivateGy    = 0;
expControl.sequence.deactivateGz    = 0;
expControl.sequence.deactivateSS    = 0; % deactivate slice selection

% motion related
expControl.motionSpecs.pattern      = 0;
expControl.motionSpecs.rotAngle     = 0;
expControl.motionSpecs.rotFreq      = 0;
expControl.motionSpecs.transMag     = 0;
expControl.motionSpecs.transFreq    = 0;
expControl.motionSpecs.transAxis    = 0;

% standard and other
expControl.advancedNotifications   = 0;
    
if fillData
    %% fill data with DB info
    try
        
        if strcmp(application,'sbr')
                       
            expControl.application                  = application;
            expControl.debug.errorFolder            = sessionData.errorFolder;            
            expControl.simulation.kernelFolder      = sessionData.kernelFolder;
            expControl.simulation.kernelPtx         = sessionData.kernelPtx;
            expControl.simulation.simulationEngine  = sessionData.simulationEngine;
            expControl.simulation.numGPU            = sessionData.numGPU;
            expControl.simulation.threads           = sessionData.threads;
            expControl.simulation.blocks            = sessionData.blocks;
            expControl.simulation.precision         = sessionData.precision;
            
            expControl.debug.debugMode              = sessionData.debugMode;
            expControl.debug.debugFile              = sessionData.debugFile;
            
            %% Experiment name
            expControl.name = sprintf('%s_%s', application, timeStamp);            

        else
            
            %% session data
            expControl.application      = sessionData.application; % 'edutool'
            expControl.versionNum       = sessionData.versionNum;
            expControl.instanceID       = sessionData.instanceID;
            expControl.courseID         = sessionData.courseID;
            expControl.connLocalDB      = sessionData.connLocalDB;          
            
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
                '''advanced_notifications'') ORDER BY FIELD ',...
                '(name,''performance_gridZ_sliceThickness'',''PDinhomogeneity'',',...
                '''isotropic_grid'',''analytical_simulation'',''fast_algorithm'',',...
                '''coil'',''isocenter'',''motion'',''motion_rot_angle'',',...
                '''motion_rot_freq'',''motion_trans_magn'',''motion_trans_freq'',',...
                '''motion_trans_axis'',''kernel_threads'',''kernel_blocks'',',...
                '''gridSizeOption'',''activate_cs'',''dwell_time'',',...
                '''coil_basic_or_advanced'',''deactivateGx'',''deactivateGy'',',...
                '''deactivateGz'',''simulationMode'',''simulationKernel'',',...
                '''advanced_notifications'')'];
            % query
            sqlQueryResults    = exec(sessionData.connLocalDB, sqlQuery);
            perfTypeInfo       = fetch(sqlQueryResults);

            %% update expControl (notice numbers of Data not in order!)

            % model related
            expControl.model.useSliceThickness  = str2double(perfTypeInfo.Data{1,1});
            expControl.model.pdInhomogeneity    = str2double(perfTypeInfo.Data{2,1});
            expControl.model.isotropicGrid      = str2double(perfTypeInfo.Data{3,1});
            expControl.model.coilType           = perfTypeInfo.Data{6,1};
            expControl.model.zIsocenter         = str2double(perfTypeInfo.Data{7,1})/100;
            if (expControl.courseID == 6)
                % MIDA, fix isocenter, need to correct in Front End
                expControl.model.zIsocenter = 0.12;
            end
            expControl.model.coilMode           = perfTypeInfo.Data{19,1};  % It refers to the option coil_basic_or_advanced
            expControl.model.gridSizeOption     = perfTypeInfo.Data{16,1};
            switch lower(expControl.model.gridSizeOption)
                case '0p5mm'
                    expControl.model.gridStep   = [0.5, 0.5, 0.5]*1e-3;
                otherwise
                    expControl.model.gridStep   = [1.0, 1.0, 1.0]*1e-3;
            end

            % simulation related
            expControl.simulation.analyticalSim     = str2double(perfTypeInfo.Data{4,1}); % full analytical
            expControl.simulation.simulationEngine  = perfTypeInfo.Data{23,1}; % numerical / analytical
            expControl.simulation.simulationKernel  = perfTypeInfo.Data{24,1}; % latest / stable

            expControl.simulation.threads           = str2double(perfTypeInfo.Data{14,1});
            expControl.simulation.blocks            = str2double(perfTypeInfo.Data{15,1});
            expControl.simulation.activateCS        = str2double(perfTypeInfo.Data{17,1});

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

        end
        
    catch ME
        ME.identifier;
        ME.message;
        
        %% inform user
        msg = ['There was a connectivity error. ',...
            'The error has been reported for further review. ',...
            'Please try again and if error persist '...
            'contact your administrator.'];
        eduTool.frontend.errorAndDBentry(sessionData.connLocalDB, msg, ...
            'cancelled-error',0,0);
        
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
    
end
