function [ok] = simulatorInstance(jobNum, nats)
%
% SERVICES.SIMULATORINSTANCE
%
%	Mock up of the SIMULATOR manager service
%
% INPUT
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'services.simulatorInstance';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

tTotal = tic();
fprintf(1, '\n STARTING SIMULATOR INSTANCE for job %d of %s',...
    jobNum, nats.experimentName);
fprintf(1, '\n');

%% Extract information for nats message
% misc data, that we do not need in this case
timeStamp       = nats.timeStamp;
userID          = nats.userID;
courseID        = nats.courseID;
experimentID    = nats.experimentID;
% DB comm is needed here, different comm channel in services
connLocalDB     = nats.connLocalDB;
redisPath       = nats.redisPath; % mock up of redis
% to load data from json file
experimentName  = nats.experimentName;

%% LOAD EXPERIMENT DATA FROM JSON FILE IN REDIS
experimentFile = sprintf('%s%s.json',redisPath,experimentName);
fid = fopen(experimentFile,'r');
experimentData = jsondecode(fread(fid,inf,'*char').');
fclose(fid);

% extract data
acquisition = experimentData.acquisition;
expControl  = experimentData.expControl;
% upgrade the expControl with the DB connectivity
expControl.connLocalDB = connLocalDB;

%% EXTRACT SIMULATION MODEL FROM REDIS:
load(sprintf('%s%s-simjob%d.mat',redisPath,experimentName,jobNum));

%% EXTRACT SEQUENCE AND MOTION FROM REDIS:
load(sprintf('%s%s-sequence.mat',redisPath,experimentName));
load(sprintf('%s%s-motion.mat',redisPath,experimentName));

%% RUN THE SIMULATION
timeSolution = [];
[timeSolution] = simulator.runSimulation( ...
    simJob.model, pulseSequence, motionModel, expControl, timeSolution );

%% SAVE SIMULATION RESULT TO REDIS
% add data to keep track of slice and part
timeSolution.sliceNum = simJob.sliceNum;
timeSolution.numParts = simJob.numParts;
timeSolution.partNum  = simJob.partNum;
% save
save(sprintf('%s%s-simresult%d.mat',redisPath,experimentName,jobNum),...
        'timeSolution', '-v7.3');

%% REPORT OK
ok = 1;
fprintf(1, '\n');
fprintf(1, '\n SIMULATOR INSTANCE for job %d on %s DONE',...
    jobNum, experimentName);
fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
fprintf(1, '\n');
