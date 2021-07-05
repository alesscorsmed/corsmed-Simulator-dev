function runStandAloneApp()
%
% EDUTOOL.RUNSTANDALONEAPP
%
%	Runs the stand alone EduTool V2 (modular design)
%
%
%========================  CORSMED AB © 2020 ==============================
%

rehash;
%% general data
APP     = 'edutool';
VERSION = 'S2V2.1.1.R20210302';
% Initiate the appSpecs structure
appSpecs.application    = APP;
appSpecs.approach       = '';
appSpecs.mode           = '';
appSpecs.version        = VERSION;
% time
tBoot   = tic();


%% initialize a pool with number of GPUs available
tParpool = tic();
CPUperGPU = 1;
edtPool = [];
gpuPool = [];
myCluster = parcluster('local');
numCPUs = myCluster.NumWorkers;
numGPUs = gpuDeviceCount;
if numGPUs > 0
    % initialize parPool if it does not exist
    numCPUs = min(numCPUs,CPUperGPU*numGPUs); % limit max of CPUs
    if isempty(gcp('nocreate'))
        edtPool = parpool(numCPUs, 'IdleTimeout', Inf);
    else
        edtPool = gcp();
    end
    % create a gpuPool with the GPU devices
    if edtPool.NumWorkers > 1
        spmd
            % each worker selects its own gpu
            gpuPool = gpuDevice(rem(labindex,numGPUs)+1);
        end
    else
        % single GPU, assign to current CPU
        gpuPool{1} = gpuDevice(labindex);
    end
else
    MException('EduTool:BadInstance', 'Instance has not available GPUs');
end
fprintf(1, '\n Parallel Pool initialized with %d Workers', edtPool.NumWorkers);
fprintf(1, '\n  Initialization Time  %.3fs', toc(tParpool));
fprintf(1, '\n');
%% initialization: AWS instance & Session attributes, connect to DB
[sessionData, tagsStruct] = eduTool.setup.initializeInstance(appSpecs);

%% enable to run
if isfield(sessionData,'connLocalDB')
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
end

%% load course info
[anatomicalModel, coilSystem] = ...
    eduTool.setup.loadCourseData( sessionData, tagsStruct );

%% BackEnd is ready for an Experiment
eduTool.frontend.updateScannerStatus(sessionData, ...
    'The virtual MR scanner is now ready...');
if isfield(sessionData,'connLocalDB')
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);
end
fprintf(1, '\n Scanner ready' );
fprintf(1, '\n  Booting Time      %.3fs', toc(tBoot));
fprintf(1, '\n\n\n eduTool ready\n');

%% Wait for an Experiment to run from FrontEnd
waitForExperiment = 1;
previousExp.coilConfirm = 0;
previousExp.tag         = '';
while waitForExperiment
    
    %% Check for a termination event
    if isfield(sessionData,'connLocalDB')
        restartMCR = eduTool.frontend.checkTerminationEvent(sessionData.connLocalDB);
        if restartMCR
            eduTool.frontend.updateScannerStatus(sessionData, ...
                'Shutting down V2 simulator...');
            break;
        end
    end
    
    %% Returns the necessary information for the experiment from FrontEnd
    [expReady,expControl,acquisition,sessionData] = ...
        eduTool.setup.fetchExperimentData(sessionData,APP);
    
    %% check if experiment is ready
    if ~expReady
        continue
    end
    
    %% start experiment >>
    
    %% new name
    experimentTAG = sprintf('U%d_C%d_E%d_P%d_UI%d',...
        expControl.userID, ...
        expControl.courseID, ...
        expControl.experimentID,...
        expControl.pulseqID,...
        acquisition.data.pulseSeqNum);
    experimentName = sprintf('%s_%s',...
        expControl.fileTimeStamp, experimentTAG);
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
    if isfield(expControl,'connLocalDB')
        eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
    end
    eduTool.frontend.updateScannerStatus(sessionData, ...
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
        % create experiment progress columns and update to 1%
        if isfield(expControl,'connLocalDB')
            exec(expControl.connLocalDB, ...
                ['INSERT INTO edt_tool_local.experiments_progress (exper_id,progress)',...
                ' VALUES (',num2str(expControl.experimentID),',0)']);
        end
        eduTool.frontend.updateExperimentProgress(expControl,'','started','')
        
        %% Update Anatomical model tissue properties from table
        [anatomicalModel] = anatomical.updateTissueProperties(...
            expControl,anatomicalModel);
        %% Update Coil Selection
        [coilSystem] = coils.updateCoilSelection(expControl,anatomicalModel,coilSystem);
        %% Verify coil parallelization and acceleration
        if (previousExp.coilConfirm) && strcmp(experimentTAG,previousExp.tag)
            % if same experiment and was confirmed, do not check
            previousExp.coilConfirm = 0;
            previousExp.tag         = '';
        else
            % need confirmation
            previousExp.coilConfirm   = 1;
            previousExp.tag           = experimentTAG;
            % if passes the check, confirmation will be set to 0
            [previousExp.coilConfirm] = coils.checkAcceleration(coilSystem, acquisition, expControl);
        end
            
        %% upgrade expControl with parPool and gpuPool
        [expControl,edtPool,gpuPool,numGPUs] = ...
            eduTool.multiGPU.setParallelization( ...
                    expControl,edtPool,gpuPool,numGPUs);
        %% call the execution
        expStats = eduTool.run.experimentExecution(...
            anatomicalModel, coilSystem, acquisition, expControl);
        
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
        %% save stats
        experimentDataFile = sprintf('%s%s.json',...
            expControl.folderSystem.statsFolder, ...
            expControl.name);
        fprintf(1, '\n');
        fprintf(1, '\n Saving statistics ... ');
        tJson = tic();
        % save experimentData to json file
        fid = fopen(experimentDataFile,'w');
        fwrite(fid, jsonencode(expStats), 'char');
        fclose(fid);
        fprintf(1, ' done -- elapsed time %.3fs', toc(tJson));
        fprintf(1, '\n  File : %s', experimentDataFile);
        fprintf(1, '\n');
        
        %% update status to 'finished'
        eduTool.frontend.updateExperimentProgress(expControl,'100','finished','')
        
        %% print result and stats
        fprintf(1, '\n');
        fprintf(1, '\n IMAGING EXPERIMENT %s DONE', experimentName);
        fprintf(1, '\n  User       ID  : %d', expControl.userID);
        fprintf(1, '\n  Course     ID  : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID  : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID  : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence #  : %d', acquisition.data.pulseSeqNum);
        fprintf(1, '\n  Instance   ID  : %s', expControl.instanceID);
        fprintf(1, '\n    # CPUs       : %d', expStats.numCPUs);
        fprintf(1, '\n    # GPUs       : %d (%s)', expStats.numGPUs, expStats.typeGPU); 
        fprintf(1, '\n  # Rx Coils     : %d', expStats.numRxCoils);
        fprintf(1, '\n  # Samples      : %d', expStats.numRxs);
        fprintf(1, '\n  # Time Steps   : %d', expStats.numSteps);
        fprintf(1, '\n  # Voxels       : %d', expStats.numVoxels);
        fprintf(1, '\n  # Slices       : %d', expStats.numSlices);
        fprintf(1, '\n  # Simulations  : %d', expStats.numJobs);
        fprintf(1, '\n  Sequence Type  : %s', expStats.seqFamily);
        fprintf(1, '\n  Sequence TIRL  : %.1fs', expStats.seqTIRL);
        fprintf(1, '\n  RUN      time  : %.1fs \t(%.1fX w.r.t TIRL)', expStats.timeTotal, expStats.seqTIRL/expStats.timeTotal);
        fprintf(1, '\n  DICOM    time  : %.1fs \t(%.1f%%)', expStats.timeDicom, 100*expStats.timeDicom/expStats.timeTotal);
        fprintf(1, '\n  RECON    time  : %.1fs \t(%.1f%%)', expStats.timeRecon, 100*expStats.timeRecon/expStats.timeTotal);
        fprintf(1, '\n  SLICER   time  : %.1fs \t(%.1f%%)', expStats.timeSlicer, 100*expStats.timeSlicer/expStats.timeTotal);
        fprintf(1, '\n  ENGINE   time  : %.1fs \t(%.1f%%)', expStats.timeEngine, 100*expStats.timeEngine/expStats.timeTotal);
        fprintf(1, '\n    Performance  : %.1f ps/voxel/step', 1e12*expStats.timeEngine/expStats.numVoxels/expStats.numSteps);
        fprintf(1, '\n  ELAPSED  time  : %.1fs',toc(tExperiment));
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
                if isfield(expControl,'connLocalDB')
                    eduTool.frontend.updateExperimentStatus(...
                        expControl.connLocalDB,expControl.experimentID,'error');
                end
                
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
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);        
        fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
        fprintf(1, '\n');
        
    end
    
    %% scanner ready to receive the data again
    eduTool.frontend.updateScannerStatus(sessionData, ...
        'The virtual MR scanner is ready');
    if isfield(sessionData,'connLocalDB')
        eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);
    end
    fprintf(1, '\n\n eduTool ready\n');
    
    %% set flag to zero to avoid inf loop
    waitForExperiment = 1;
    
end

