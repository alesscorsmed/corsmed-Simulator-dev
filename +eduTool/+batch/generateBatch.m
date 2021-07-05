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

%% clear and prepare
clear all; close all; clc;

%% prepare data to connect with Front End
rehash
VERSION = 'V2_9_19';
disp(['VERSION: ',VERSION])
%% general data
APP             = 'edutool';
APPROACH        = 'standalone';
MODE            = 'test-gen';
jsontestfile    = '';
testPauseTime   = '';
simulatorSpecs  = '';
%% for test
testFolder     = '/efs-mount-point/S20/TESTS/edutool/regression'; % main folder
testTag        = 'default'; % sub-folder for different experiment configurations
testNum        = 0; % number of the test

%% test control struct
FAIL = 0;
if exist(testFolder, 'dir')
    testControl.testFolder = testFolder;
    % create subfolders folders
    if mkdir(testFolder,'inputs')
        testControl.inputsFolder = sprintf('%s/inputs',testControl.testFolder);
    else
        FAIL = 1;
    end
    if mkdir(testFolder,'golden')
        testControl.goldenFolder = sprintf('%s/golden',testControl.testFolder);
    else
        FAIL = 1;
    end
    if mkdir(testFolder,'result')
        testControl.resultFolder = sprintf('%s/result',testControl.testFolder);
    else
        FAIL = 1;
    end
else
    FAIL = 1;
end
if FAIL
    msg = sprintf('Batch folder %s does not exist or cannot create subfolders',...
        testFolder);
    ME = MException('error:generateRegressionTests', '%s', msg);
    throw(ME);
end
testControl.caseTag = testTag;
testControl.version = VERSION;

%% initialize instance and app data
appSpecs = eduTool.setup.initializeApplication(APP,APPROACH,MODE,VERSION,...
    jsontestfile,testPauseTime);

tBoot   = tic();

%% initialization: AWS instance & Session attributes, connect to DB
[sessionData, tagsStruct] = eduTool.setup.initializeInstance(appSpecs);

%% make additions in the sessionData structure
sessionData.jsontestfile    = appSpecs.jsontestfile;
sessionData.testPauseTime   = str2num(appSpecs.testPauseTime);

%% enable to run
if isfield(sessionData,'connLocalDB')
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
end

%% load course info
[anatomicalModel, coilSystem, sessionData] = ...
    eduTool.setup.loadCourseData( sessionData, tagsStruct );

%% BackEnd is ready for an Experiment
eduTool.frontend.updateScannerStatus(sessionData, ...
    'The virtual MR scanner is now ready...');
if isfield(sessionData,'connLocalDB')
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);
end
fprintf(1, '\n Scanner ready' );

%% Wait for an Experiment to run from FrontEnd
waitForExperiment = 1;
while waitForExperiment
    
    %% Returns the necessary information for the experiment from FrontEnd
    [expReady,expControl,acquisition,sessionData] = ...
        eduTool.setup.fetchExperimentData(sessionData,APP,APPROACH,MODE);
    
    %% check if experiment is ready
    if ~expReady
        continue
    end
    
    %% Experiment
    experimentTAG = sprintf('U%d_C%d_E%d_P%d_UI%d',...
        expControl.userID, ...
        expControl.courseID, ...
        expControl.experimentID,...
        expControl.pulseqID,...
        acquisition.data.pulseSeqNum);
    experimentName = sprintf('%s_%s',...
        expControl.fileTimeStamp, experimentTAG);
    
    tExperiment = tic();
    fprintf(1, '\n STARTING IMAGING EXPERIMENT %s', experimentName);
    fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
    fprintf(1, '\n  User       ID : %d', expControl.userID);
    fprintf(1, '\n  Course     ID : %d', expControl.courseID);
    fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
    fprintf(1, '\n  Pulse Seq. ID : %d', expControl.pulseqID);
    fprintf(1, '\n  UI Sequence # : %d', acquisition.data.pulseSeqNum);
    fprintf(1, '\n');
    expControl.name = experimentName;
    
    %% test name
    testName   = eduTool.batch.generateExperimentName(...
        acquisition.data,expControl);
    % add version, case and time stamp
    experimentName = sprintf('%s_%s_%s_%s',...
        testControl.version, testControl.caseTag, ...
        batchName, expControl.fileTimeStamp);
    expControl.name = experimentName;

    %% encode experiment in JSON format
    % data into structure
    experimentData = [];
    experimentData.expControl  = expControl;
    experimentData.acquisition = acquisition;
    % do not save connLocalDB
    experimentData.expControl.connLocalDB   = [];
    experimentData.expControl.application   = 'test';
    % convert and saveto json
    experimentData = jsonencode(experimentData);
    
    %% disable run button and start experiment
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
    eduTool.frontend.updateScannerStatus(sessionData, ...
        'Saving experiment to batch...');
    
    %% report experiment done
    tExperiment = toc(tExperiment);
    
    %% save experiment 
    experimentDataFile = sprintf('%s/%s.json',...
        testControl.inputsFolder, ...
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
    eduTool.frontend.updateExperimentStatus(...
        expControl.connLocalDB,expControl.experimentID,'finished');
    
    %% print result
    fprintf(1, '\n');
    fprintf(1, '\n DONE');
    fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
    fprintf(1, '\n');

    %% scanner ready to receive the data again
    eduTool.frontend.updateScannerStatus(sessionData, ...
        'The virtual MR scanner is ready');
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);
    
    %% set flag to zero to avoid inf loop
    waitForExperiment = 1;
    
end

