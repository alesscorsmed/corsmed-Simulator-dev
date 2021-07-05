sessionData.userID = user_id;
sessionData.userGroupID = usergroup_id;

sessionData.instanceID = instance_id;
sessionData.courseID = course_id;
sessionData.backendVersion = versionNum;

sessionData.experimentID = experiment_id;
sessionData.experimentInfo = exper_info;

% DB
sessionData.connLocalDB = conn_localdb;
sessionData.localDB = localdb;
sessionData.connRemoteDB = conn_remotedb;
sessionData.remoteDBID = remoteDB_id;
sessionData.talkToRemoteDB =talkToRemoteDB;

% folders
sessionData.kernelFolder = kernelFolder;
sessionData.commonFolder =commonFolder;
sessionData.resultsFolder =resultsFolder;
sessionData.gadgetronFolder =gadgetronFolder;
sessionData.gadgetronResultsFolder =gadgetronResultsFolder;
sessionData.gadgetronISMRMRDFolder =gadgetronISMRMRDFolder;
sessionData.analyticReconFolder =analyticReconFolder;

% mode
sessionData.testMode =testMode;

% sequence data
sessionData.pathPulseSeq =pathPulseSeq;
sessionData.pulseqID =pulseq_id;

% recon info
sessionData.reconstructor =reconstructor;
sessionData.reconInfo =recon_info;

sessionData.outerFOVratio =outerFOVratio;
sessionData.stackStruct =stack_struct;
sessionData.defaultGridStep =default_gridStep;

% tags
sessionData.tagParfeval = tag_parfeval;
sessionData.tagStruct = tags_struct;

% initialize MR system
mrSystem = data.mrSystem.initialize();
% dummys for now
anatomicalModel = [];
coilSystem      = [];

% call the runExperiment
[expInfo] = modular.singleGPU.runExperiment(...
    sessionData, anatomicalModel, mrSystem, coilSystem);



