function [nats] = simulatorManager(nats)
%
% SERVICES.SIMULATORMANAGER
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
functionName = 'services.simulatorManager';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

tTotal = tic();
fprintf(1, '\n STARTING SIMULATOR MANAGER on %s', nats.experimentName);
fprintf(1, '\n');

%% ALL INSTANCE SCALING AND CREATION
% here just execute the jobs added to the queue sequentially
jobsDone = 0;
for jobNum = 1:nats.numJobs
    % call the instance to simulate
    [ok] = services.simulatorInstance(jobNum, nats);
    jobsDone = jobsDone + ok;
end

%% REPORT OK
nats.simulatorService = jobsDone;
fprintf(1, '\n');
fprintf(1, '\n SIMULATOR MANAGER on %s DONE', nats.experimentName);
fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
fprintf(1, '\n');
