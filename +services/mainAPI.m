% 
%  main API entry point (mockup)
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

%% initialization: AWS instance & Session attributes, connect to DB
[sessionData, tagsStruct] = eduTool.setup.initializeInstance();

% run button not ready yet
eduTool.frontend.updateScannerStatus(sessionData, ...
    'The virtual MR scanner is booting ...');
eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);

%% TO REMOVE: only for testing purposes
% change folders in sessionData
sessionData.folderSystem.anatomicalModelFolder  = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/anatomical_models/';
sessionData.folderSystem.coilModelFolder        = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/coil_models/';
sessionData.folderSystem.errorFolder            = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/errors/';
sessionData.folderSystem.experimentFolder       = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/experiments/';

%% MOCK UP OF NATS
% DB conn is needed here, different comm channel in services
nats.connLocalDB = sessionData.connLocalDB;
nats.redisPath = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/REDIS/'; % mock up of redis
nats.s3Path    = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/S3/'; % mock up of S3

%% BackEnd is ready for an Experiment
eduTool.frontend.updateScannerStatus(sessionData, ...
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
    
    %% disable run button and start experiment
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
    eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
        'Running experiment...');
    
    %% Returns the necessary information for the experiment from FrontEnd
    [expControl, acquisition] = eduTool.setup.fetchExperimentData( ...
        experiment, sessionData );
    
    % upgrade with some temporaty fields
    expControl.useOldSequence = 0;
        
    %% change status to 'started'
    eduTool.frontend.updateExperimentStatus(expControl.connLocalDB,...
        expControl.experimentID,'started')
    
    %% create experiment progress columns and update to 1%
    exec(expControl.connLocalDB, ...
        ['INSERT INTO edt_tool_local.experiments_progress (exper_id,progress)',...
        ' VALUES (',num2str(expControl.experimentID),',1)']);
    
    %% Experiment name
    experimentName = sprintf('%s-U%d-C%d-E%d',...
        expControl.timeStamp, ...
        expControl.userID, ...
        expControl.courseID, ...
        expControl.experimentID);
    
    %% encode experiment in JSON format
    % data into structure
    experimentData = [];
    experimentData.expControl  = expControl;
    experimentData.acquisition = acquisition;
    % do not save connLocalDB
    experimentData.expControl.connLocalDB   = [];
    % convert and save to json
    experimentData = jsonencode(experimentData);
    experimentDataFile = sprintf('%s%s.json',...
        nats.redisPath, experimentName);
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
    
    %% ADD EXPERIMENT TO NATS
    nats.timeStamp      = expControl.timeStamp;
    nats.userID         = expControl.userID;
    nats.courseID       = expControl.courseID;
	nats.experimentID   = expControl.experimentID;
    nats.experimentName = experimentName;
    
    try
        
        %% Start the Experiment
        tExperiment = tic();
        fprintf(1, '\n STARTING IMAGING EXPERIMENT %s', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n');
        
        %% SLICER SERVICE
        [nats] = services.slicerManager(nats);
        
        %% SIMULATION SERVICE
        [nats] = services.simulatorManager(nats);
        
%         %% RECONSTRUCTION SERVICE
%         [nats] = services.reconstructionManager(nats);

        %% report experiment done
        tExperiment = toc(tExperiment);
        
        %% save experiment as success to S3
        experimentDataFile = sprintf('%s%s.json',...
            nats.s3Path, experimentName);
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
        
        %% print result
        fprintf(1, '\n');
        fprintf(1, '\n IMAGING EXPERIMENT %s DONE', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
    
    catch ME
        
        %% catch the error
        ME.identifier;
        ME.message;
        
        %% if cancelled by user, no backend report
        if strcmp(ME.message,'CANCELLED-BY-USER')
            
            %% save experiment as cancelled in experiment
            experimentDataFile = sprintf('%s%s-CANCELLED.json',...
                expControl.folderSystem.errorFolder, experimentName);
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
            
        else
            
            %% update status to 'Error' if not due to user cancellation
            eduTool.frontend.updateExperimentStatus(...
                expControl.connLocalDB,expControl.experimentID,'error');

            %% save experiment as error for debugging
            experimentDataFile = sprintf('%s%s-ERROR.json',...
                sessionData.folderSystem.errorFolder, experimentName);
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
                eduTool.frontend.updateExperimentProgress(expControl,'','cancelled-error',msg);
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

