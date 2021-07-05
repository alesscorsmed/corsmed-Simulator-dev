function [] = initializeEduTool()
%
% EDUTOOL.INITIALIZEEDUTOOL
%
%	initializes the EduTool and waits for experiment
%
% INPUT
%
%   none
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%


%% initialize AWS instance & Session attributes
[instanceData] = data.instance.initialize();
[instanceData,tagsStruct] = aws.ec2.readInstanceTags();


%% connect to DB
[expControl] = data.expControl.initialize();
[sessionData] = data.sessionData.initialize();
[expControl,sessionData] = backend.connectDatabase(tagsStruct,instanceData,expControl);


%% Load the human anatomical model 
%initialize Anatomical Model
%[anatomicalModel] = data.anatomicalModel.initialize();

%first load anatomical model data from session
%sessionData.AnatModel
sessionData = eduTool.findAnatomicalModelPath(sessionData,instanceData.Courseid);


%% Update DB that backend has started
eduTool.frontend.updateBackendStartCounts(expControl.connLocalDB)


%% Parpool Connection if necessary
if(str2num(instanceData.parfeval) == 1)
    parpoolConn = backend.parpoolConnection(expControl.connLocalDB,expControl.localDB);
else
    parpoolConn = 0;
end

%% Load Coil Profiles

%% Unitaty MR-safety and SAR pre-computation

%% PD inhomogeneity

%% BackEnd is ready for an Experiment
eduTool.frontend.updateScannerStatus(expControl.connLocalDB,'The virtual MR scanner is now ready...');
eduTool.frontend.runButtonEnabled(expControl.connLocalDB,1);

%% Wait for an Experiment to run from FrontEnd
    while 1

        
        experimentReady = eduTool.frontend.expectExperiment(expControl.connLocalDB);      
        if strcmp(experimentReady.Data,'No Data')
            continue
        end
        
        
        %% update Coil and MR safety parameters if needed
            %[coils_struct,MRsafety_struct] = sar.updateMRSafety(conn_localdb,...
            %      anatHumanModel_struct,coils_struct,MRsafety_struct,course_id);
            %  edt_db_updateScannerStatus(conn_localdb,'The virtual MR scanner is ready.');       
        clc
        close all

        
        %% Returns the necessary information for the experiment from Front
        [expControl] = eduTool.frontend.fetchExperimentInfo(experimentReady,expControl,sessionData);             
        
        
        %% Start the Experiment
        disp('THIS IS AN IMAGING EXPERIMENT')
        info4user = backend.runExperiment(sessionData,experimentInfo,expControl,tagsStruct,parpoolConn);
        
        
        
    end