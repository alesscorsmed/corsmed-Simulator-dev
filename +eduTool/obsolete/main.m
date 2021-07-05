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

%% initialization: AWS instance & Session attributes, connect to DB
[sessionData, tagsStruct] = eduTool.setup.initializeInstance();

eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);


%% Create Parallel Pool
p = eduTool.setup.openParpool();

%% TO REMOVE: only for testing purposes
% change folders in sessionData
sessionData.folderSystem.anatomicalModelFolder  = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/anatomical_models/';
sessionData.folderSystem.coilModelFolder        = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/coil_models/';
sessionData.folderSystem.errorFolder            = '/efs-mount-point-MATLAB/EduTool-Jorge/FILES/MODULAR/errors/';


%% load course info
[anatomicalModel, mrSystem, coilSystem] = ...
    eduTool.setup.loadCourseData( sessionData );

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
    
    % disable run button and start experiment
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,0);
    eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
    'Running experiment...');
        
    %% Returns the necessary information for the experiment from FrontEnd
    [expControl, acquisition] = eduTool.setup.fetchExperimentData( ...
        experiment, sessionData );
    
    %% compute Noise levels
    [acquisition] = noise.calculateNoiseLevel( acquisition, ...
        expControl, mrSystem );
    
    %% Update Coil Selection
    [coilSystem] = coils.updateCoilSelection(...
        expControl, anatomicalModel, coilSystem );
    
    %% Start the Experiment
    tExperiment = tic();
    fprintf(1, '\n STARTING IMAGING EXPERIMENT');
    fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
    fprintf(1, '\n  User       ID : %d', expControl.userID);
    fprintf(1, '\n  Course     ID : %d', expControl.courseID);
    fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
    fprintf(1, '\n');
    
    %% select simulation engine
    expControl.simulation.simulationEngine = 'analytical'; % analytical / phasor / diffusion
    
    %% in case we want to use old sequences, only analytical works
    expControl.useOldSequence = 0;
    if expControl.useOldSequence
        expControl.simulation.simulationEngine = 'analytical';
    end
    
    % optimization flags (only for new sequences)
    expControl.sequence.timeCompression = 1; % time compression
    expControl.sequence.minContextExc   = 1; % reduces number of kernel calls
    
    % shows back-end progress bar
    %  CAUTION: impacts performance
    %           to time performance:
    %               set expControl.debug.waitBarBE = 0;
    %               set expControl.debug.debugMode = 0;
    expControl.debug.waitBarBE = 0;
    expControl.debug.debugMode = 1;
    
    % encoding reorder
    acquisition.data.encOrder = 'std'; % std / centric / edge
    
    %% call the execution
    info4user = eduTool.run.experimentExecution(...
        anatomicalModel, coilSystem, mrSystem, ...
        acquisition, expControl, tagsStruct, parpoolConn, sessionData);
    
    %% report experiment done
    tExperiment = toc(tExperiment);
    fprintf(1, '\n');
    fprintf(1, '\n IMAGING EXPERIMENT DONE');
    fprintf(1, '\n  Instance   ID : %s', expControl.instanceID);
    fprintf(1, '\n  User       ID : %d', expControl.userID);
    fprintf(1, '\n  Course     ID : %d', expControl.courseID);
    fprintf(1, '\n  Experiment ID : %d', expControl.experimentID);
    fprintf(1, '\n  Elapsed time  : %.3fs', tExperiment);
    fprintf(1, '\n');
    
    %% scanner ready to receive the data again
    eduTool.frontend.updateScannerStatus(sessionData.connLocalDB, ...
    'The virtual MR scanner is ready');
    eduTool.frontend.runButtonEnabled(sessionData.connLocalDB,1);
    
    %% set flag to zero to avoid inf loop
    waitForExperiment = 1;
    
end