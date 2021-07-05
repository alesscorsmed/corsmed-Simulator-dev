%function [] = main()
%
% EDUTOOL.MAIN
%
%	Runs the EduTool App
%
% INPUT
%
%   none
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

%% general data
APP     = 'edutool';
VERSION = 'v20210118S2';

%% initialization: AWS instance & Session attributes, connect to DB
[sessionData, tagsStruct] = eduTool.setup.initializeInstance();
% assign application and version
sessionData.application	= APP;
sessionData.versionNum  = VERSION;

eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);

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

% %% TO REMOVE: only for testing purposes
% % change folders in sessionData
% sessionData.folderSystem.anatomicalModelFolder  = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/anatomical_models/';
% sessionData.folderSystem.coilModelFolder        = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/coil_models/';
% sessionData.folderSystem.errorFolder            = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/errors/';
% sessionData.folderSystem.experimentFolder       = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/experiments/';

%% load course info
[anatomicalModel, mrSystem, coilSystem] = ...
    eduTool.setup.loadCourseData(sessionData,tagsStruct);

%% Parpool Connection if necessary
if(sessionData.parfeval == 1)
    parpoolConn = eduTool.setup.parpoolConnection(...
        sessionData.connLocalDB,sessionData.localDB);
else
    parpoolConn = 0;
end

%% BackEnd is ready for an Experiment
eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
    'The virtual MR scanner is now ready...');
eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);

%% Wait for an Experiment to run from FrontEnd
waitForExperiment = 1;
while waitForExperiment
    
    %% check for experiment
    experiment = eduTool.frontend.expectExperiment(sessionData.connLocalDB);
    if strcmp(experiment.Data,'No Data')
        continue
    end
    
    %% Returns the necessary information for the experiment from FrontEnd
    [expControl, acquisition] = eduTool.setup.fetchExperimentData( ...
        experiment, sessionData );
    
    %% compute Noise levels
    [acquisition] = noise.calculateNoiseLevel( acquisition, ...
        expControl, mrSystem );
    
    %% new name
    experimentName = sprintf('%s_U%d_C%d_E%d_P%d_UI%d',...
        expControl.fileTimeStamp, ...
        expControl.userID, ...
        expControl.courseID, ...
        expControl.experimentID,...
        expControl.pulseqID,...
        acquisition.data.pulseSeqNum);
    expControl.name = experimentName;

    %% encode experiment in JSON format
    % data into structure
    experimentData = [];
    experimentData.expControl  = expControl;
    experimentData.acquisition = acquisition;
    % do not save connLocalDB
    experimentData.expControl.connLocalDB   = [];
    % convert and saveto json
    experimentData = jsonencode(experimentData);
    
    %% disable run button and start experiment
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
    eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
        'Running experiment...');
    
    try
        
        %% Start the Experiment
        tExperiment = tic();
        fprintf(1, '\n STARTING IMAGING EXPERIMENT %s', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);        
        fprintf(1, '\n');
        
        %% change status to 'started'
        eduTool.frontend.updateExperimentStatus( ...
            expControl.connLocalDB,expControl.experimentID,'started' );
        
        %% create experiment progress columns and update to 1%
        exec(expControl.connLocalDB, ...
            ['INSERT INTO edt_tool_local.experiments_progress (exper_id,progress)',...
            ' VALUES (',num2str(expControl.experimentID),',0)']);
                
        %% Update Anatomical model tissue properties from table
        [anatomicalModel] = anatomical.updateTissueProperties(...
            expControl, anatomicalModel );
        
        %% Update Coil Selection
        [coilSystem] = coils.updateCoilSelection(...
            expControl, anatomicalModel, coilSystem );
        
        %% select simulation engine
        %expControl.simulation.simulationEngine = 'analytical'; % analytical / numerical / phasor / diffusion
        
        %% in case we want to use old sequences, only analytical works
        expControl.useOldSequence = 0;
        if expControl.useOldSequence
            expControl.simulation.simulationEngine = 'analytical';
        end
        
        % shows back-end progress bar
        %  CAUTION: impacts performance
        %           to time performance:
        %               set expControl.debug.waitBarBE = 0;
        %               set expControl.debug.debugMode = 0;
        expControl.debug.waitBarBE = 0;
        expControl.debug.debugMode = 1;
        
        %% call the execution
        info4user = eduTool.run.experimentExecutionV2(...
            anatomicalModel, coilSystem, mrSystem, ...
            acquisition, expControl, tagsStruct, parpoolConn);
        
        %% report experiment done
        tExperiment = toc(tExperiment);
        
        %% save experiment as success
        experimentDataFile = sprintf('%s%s.json',...
            expControl.folderSystem.experimentFolder, ...
            expControl.name);
        fprintf(1, '\n');
        fprintf(1, '\n Saving experiment ... ');
        tJson = tic();
        % save experimentData to json file
        fid = fopen(experimentDataFile,'w');
        fwrite(fid, experimentData, 'char');
        fclose(fid);
        fprintf(1, ' done -- elapsed time %.3fs', toc(tJson));
        fprintf(1, '\n  File : %s', experimentDataFile);
        fprintf(1, '\n');
        
        %% update status to 'finished'
        eduTool.frontend.updateExperimentStatus( ...
            expControl.connLocalDB,expControl.experimentID,'finished');
        
        %% print result
        fprintf(1, '\n');
        fprintf(1, '\n IMAGING EXPERIMENT %s DONE', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);        
        fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
        
    catch ME
        
        %% catch the error
        ME.identifier;
        ME.message;
        
        switch ME.message
            
            case 'CANCELLED-BY-USER' %% if cancelled by user, no backend report
                
                %% save experiment as cancelled in experiment
                experimentDataFile = sprintf('%s%s_CANCELLED.json',...
                    expControl.folderSystem.errorFolder, ...
                    expControl.name);
                % save experimentData to json file
                fprintf(1, '\n');
                fprintf(1, '\n Saving CANCELLED experiment ... ');
                tJson = tic();
                % save experimentData to json file
                fid = fopen(experimentDataFile,'w');
                fwrite(fid, experimentData, 'char');
                fclose(fid);
                fprintf(1, ' done -- elapsed time %.3fs', toc(tJson));
                fprintf(1, '\n  File : %s', experimentDataFile);
                fprintf(1, '\n');
                
                %% error message
                errorMessage = sprintf(['%s - User (%d) - Exper (%d)',...
                    '\n CANCELLED-BY-USER',...
                    '\n Data for replication saved in %s'], ...
                    expControl.timeStamp, ...
                    expControl.userID,...
                    expControl.experimentID,...
                    experimentDataFile);
                
            case 'MESSAGE-BY-PLATFORM' %% platform error
                
                %% error message
                errorMessage = sprintf(['%s - User (%d) - Exper (%d)',...
                    '\n MESSAGE-BY-PLATFORM: revise experiment parameters'], ...
                    expControl.timeStamp, ...
                    expControl.userID,...
                    expControl.experimentID);
                
            otherwise %% actual backend error
                
                %% update status to 'Error' if not due to user cancellation
                eduTool.frontend.updateExperimentStatus(...
                    expControl.connLocalDB,expControl.experimentID,'error');
                
                %% save experiment as error for debugging
                experimentDataFile = sprintf('%s%s_ERROR.json',...
                    expControl.folderSystem.errorFolder, ...
                    expControl.name);
                % save experimentData to json file
                fprintf(1, '\n');
                fprintf(1, '\n Saving ERROR experiment ... ');
                tJson = tic();
                % save experimentData to json file
                fid = fopen(experimentDataFile,'w');
                fwrite(fid, experimentData, 'char');
                fclose(fid);
                fprintf(1, ' done -- elapsed time %.3fs', toc(tJson));
                fprintf(1, '\n  File : %s', experimentDataFile);
                fprintf(1, '\n');
                
                %% send error in connection to DB for backend
                errorMessage = sprintf(['%s - User (%d) - Exper (%d)',...
                    '\n Error in function %s() at line %d.',...
                    '\n Error Message: %s',...
                    '\n Data for error replication saved in %s'], ...
                    expControl.timeStamp, ...
                    expControl.userID,...
                    expControl.experimentID,...
                    ME.stack(1).name,ME.stack(1).line,ME.message,...
                    experimentDataFile);
                
                
                %% notify in slack if not developers
                if ~ismember(expControl.userID,[790,933,1139])
                    
                    %% slack message
                    eduTool.frontend.notifyAdminForErrors(expControl.connLocalDB,...
                        errorMessage,expControl.instanceID)
                    
                    %% inform user
                    msg = ['The experiment stopped and quit unexpectedly. ',...
                        'This issue has been reported to CORSMED for further review. ',...
                        'If you were so kind as to help us in this process, ',...
                        'please send us feedback (top right tab) with details about the experiment. ',...
                        'We apologize and will try to address this issue in the coming release.'];
                    eduTool.frontend.errorAndDBentry(expControl.connLocalDB, msg, ...
                        'cancelled-error',expControl.experimentID,expControl.pulseqID);
                end
                
        end
        
        %% display error in cmd line
        tExperiment = toc(tExperiment);
        fprintf(1, '\n');
        fprintf(1, '\n ERROR: %s', errorMessage);
        fprintf(1, '\n');
        fprintf(1, '\n IMAGING EXPERIMENT %s FAILED', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);        
        fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
        
    end
    
    %% scanner ready to receive the data again
    eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
        'The virtual MR scanner is ready');
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);
    
    %% set flag to zero to avoid inf loop
    waitForExperiment = 1;
    
end

