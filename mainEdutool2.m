function mainEdutool2(appName,jsontestfile,testPauseTime,simulatorSpecs)
%
% EDUTOOL.MAIN
%
%	Runs the EduTool App v2 (modular design) under the new architecture and
%	UI
%
% INPUT
%
%   appName: it consists of three substrings separated with the '-'
%     character. The first substring defines the APPLICATION, the second
%     substring defines the APPROACH and the third substring defines the
%     MODE. The appName can be empty and the substrings can be omitted.
%     Default case: edutool-standalone
%
%   Options:
%                                APP   -   APPROACH     |     DESCRIPTION
%                                                       |
%   1. edutool-jsonstandalone: edutool   jsonstandalone | New arch, standalone (Tested, it is working)
%   2. edutool-standalone:     edutool     standalone   | Old Arch, standalone (DEFAULT, It is working)
%   3. edutool-distributed:    edutool    distributed   | New arch, distributed 
%   4. edutool-centralized:    edutool    centralized   | New arch, centralized (Tested, It is working)
%   5. edutool-expjson:        edutool    expjson       | Read from json file
%   6. sbr:                      sbr                    | Simulation-based recon
%
%   jsontestfile: it will be used only when mode=test. It keeps the path of
%     the json file where the test experiment is stored. Currently available
%     for the edutool-jsonstandalone-test
%   testPauseTime: it will be used only when mode=test. It defines the
%     range of pause time before the simulation starts (for stress test
%     purposes). Currently available for the edutool-jsonstandalone-test
%   simulatorSpecs: available for the edutool-centralized. It defines the
%     userID, the courseID and the anatomicalID in the following form:
%     'userID-courseID-anatomicalID', for example '790-6-6'
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

rehash
vFile = dir('version*.txt');
if size(vFile,1)==1
    VERSION = vFile(1).name(9:end-4);
else
    VERSION = 'S2v20210422';
end
disp(['VERSION: ',VERSION])
%% general data
if nargin<1
    APP             = 'edutool';
    APPROACH        = 'standalone';
    MODE            = '';
    jsontestfile    = '';
    testPauseTime   = '';
    simulatorSpecs  = '';
else
    appInputs	= strsplit(appName,'-');
    APP         = appInputs{1,1};
    if size(appInputs,2)==1        
        APPROACH 	= '';
        MODE        = '';
    elseif size(appInputs,2)==2
        APPROACH 	= appInputs{1,2};
        MODE        = '';
    else        
        APPROACH 	= appInputs{1,2};
        MODE        = appInputs{1,3};
    end
    if nargin==1
        jsontestfile    = '';
        testPauseTime   = '';
        simulatorSpecs  = '';
    elseif nargin==2
        testPauseTime   = '';
        simulatorSpecs  = '';
    elseif nargin==3
        simulatorSpecs  = '';
    end
end

appSpecs = eduTool.setup.initializeApplication(APP,APPROACH,MODE,VERSION,...
    jsontestfile,testPauseTime);

tBoot   = tic();

%% initialize a pool with number of GPUs available
tParpool = tic();
CPUperGPU = 2;
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
[anatomicalModel, coilSystem, sessionData] = ...
    eduTool.setup.loadCourseData( sessionData, tagsStruct );

%% make additions in the sessionData structure
sessionData.jsontestfile    = appSpecs.jsontestfile;
sessionData.testPauseTime   = str2num(appSpecs.testPauseTime);

%% create a unique name for the matlab log files
% time stamp
timeStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
% change time stamp into file format yyyymmddHHMMSS
timeStamp = strrep(timeStamp,'-','');
timeStamp = strrep(timeStamp,':','');
timeStamp = strrep(timeStamp,' ','');
tempLogFile = strcat(timeStamp,'_',char(64+randi(84-65,1,20)),'.txt');

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
% info for confirmed experiments (to repeat)
previousExp.coilConfirm = 0;
previousExp.tag         = '';
while waitForExperiment
    try
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
            eduTool.setup.fetchExperimentData(sessionData,APP,APPROACH,MODE);

        %% check if experiment is ready
        if ~expReady
            continue
        end
        
        % Start logging
        diary(tempLogFile)

        %% Start the Experiment
        tExperiment = tic();
        
        experimentTAG = sprintf('U%d_C%d_E%d_P%d_UI%d',...
            expControl.userID, ...
            expControl.courseID, ...
            expControl.experimentID,...
            expControl.pulseqID,...
            acquisition.data.pulseSeqNum);
        experimentName = sprintf('%s_%s',...
            expControl.fileTimeStamp, experimentTAG);

        fprintf(1, '\n STARTING IMAGING EXPERIMENT %s', experimentName);
        fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
        fprintf(1, '\n  User       ID : %d', expControl.userID);
        fprintf(1, '\n  Course     ID : %d', expControl.courseID);
        fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
        fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
        fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);        
        fprintf(1, '\n');
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
        
        %% change status to 'started'
        % create experiment progress columns and update to 1%
        if isfield(expControl,'connLocalDB')
            exec(expControl.connLocalDB, ...
                ['INSERT INTO edt_tool_local.experiments_progress (exper_id,progress)',...
                ' VALUES (',num2str(expControl.experimentID),',0)']);
        end
        eduTool.frontend.updateExperimentProgress(expControl,'','started','')
        
        %% Update Anatomical model tissue properties from table
        [anatomicalModelSimul,expControl] = ...
            anatomical.updateTissueProperties(...
            expControl,anatomicalModel,acquisition.data.gamma,...
            acquisition.data.pulseSeqFamilyName);
        %% Update Coil Selection
        [coilSystem] = coils.updateCoilSelection(expControl,anatomicalModelSimul,coilSystem);
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
            anatomicalModelSimul, coilSystem, acquisition, expControl, ...
            sessionData);
        
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
        experimentDataFile = sprintf('%s%s_STATS.json',...
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
        
        %% update status to 'finished' (check status before)
        eduTool.frontend.checkExperimentStatus(expControl,expControl.experimentID)
        eduTool.frontend.updateExperimentProgress(expControl,'','finished','')
        
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
                
            otherwise %% actual backend error. In this case, notify admins
                
                msg = ['The experiment stopped and quit unexpectedly. ',...
                    'This issue has been reported to CORSMED for further review. ',...
                    'If you were so kind as to help us in this process, ',...
                    'please send us feedback (top right tab) with details about the experiment. ',...
                    'We apologize and will try to address this issue in the coming release.'];
                
                if exist('expControl','var')~=1 || isempty(expControl)
                    
                    % Notify the user regarding the current situation
                    eduTool.frontend.updateScannerStatus(sessionData,msg,'error')
                    
                else
                    %% save experiment as error for debugging
                    experimentDataFile = sprintf('%s%s_ERROR.json',...
                        expControl.folderSystem.experimentFolder, ...
                        expControl.name);  
                    
                    %% send error to admins
                    errorMessage = sprintf(['%s - User (%d) - Exper (%d)',...
                        '\n Error in function %s() at line %d.',...
                        '\n Error Message: %s',...
                        '\n Data for error replication saved in %s'], ...
                        expControl.timeStamp, ...
                        expControl.userID,...
                        expControl.experimentID,...
                        ME.stack(1).name,ME.stack(1).line,ME.message,...
                        experimentDataFile); 
                    
                    %%
                    if isfield(expControl,'connLocalDB')
                        eduTool.frontend.updateExperimentStatus(...
                            expControl.connLocalDB,expControl.experimentID,'error');
                    else
                        % Remove any special character from the messages
                        % before you send them to the users and admins
                        errorAdmins = strrep(errorMessage,'\n','');
                        errorAdmins = strrep(errorAdmins,'\r','');
                        errorAdmins = strrep(errorAdmins,'char(10)','');
                        eduTool.frontend.updateExperimentProgress(expControl,...
                            '','error',msg,errorAdmins);
                    end        

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
                    
                    %% notify in slack if not developers
                    if ~ismember(expControl.userID,[790,933,1139])                    
                        %% slack message
                        if isfield(expControl,'connLocalDB')                    
                            eduTool.frontend.notifyAdminForErrors(expControl.connLocalDB,...
                                errorMessage,expControl.instanceID)
                        end
                    end
                    
                end
                
                if isfield(sessionData,'redis')
                    % If an unexpected error occurs, empty the 
                    % EXPERIMENT_{InstanceID} redis key
                    experimentArray(1).experiment_id    = '';
                    experimentArray(1).status           = '';
                    tools.redis.redisSetJsonWrapper(sessionData.redis.R,...
                        sessionData.redis.keys.experimentsRedisKey,...
                        experimentArray);
                end
                
        end
        
        %% display error in cmd line
        if exist('expControl','var')==1 && ~isempty(expControl)            
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
        
    end
    
    eduTool.frontend.updateScannerStatus(sessionData, ...
        'The virtual MR scanner is now ready...');
        
    %% scanner ready to receive the data again
    if strcmp(APPROACH,'standalone') && isfield(sessionData,'connLocalDB')
        eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);        
    else
        %% save matlab log        
        if exist('expControl','var')~=1 || isempty(expControl)  % unexpected error            
            % time stamp
            timeStampLog = datestr(now,'yyyy-mm-dd HH:MM:SS');
            % change time stamp into file format yyyymmddHHMMSS
            timeStampLog = strrep(timeStampLog,'-','');
            timeStampLog = strrep(timeStampLog,':','');
            timeStampLog = strrep(timeStampLog,' ','');
            logDataFile = sprintf(...
                '/home/ubuntu/%s_unexpectedError_U%d_C%d_LOG.txt',...
                timeStampLog, ...
                sessionData.userID, ...
                sessionData.courseID);            
        else
            % if you make any change here, change also the path of the log
            % file stored in the jsonStructure in edutool.run.dicom
            logDataFile = sprintf('%s%s_LOG.txt',...
                expControl.folderSystem.logFolder, ...
                expControl.name);
        end
        
        % Stop logging
        diary off
        % Move the log to the correct path
        movefile(tempLogFile,logDataFile)
        
    end
    fprintf(1, '\n\n eduTool ready\n');
    
    
    %% set flag to zero to avoid inf loop
    waitForExperiment = 1;
    
end

