function [nats] = reconstructorManager(nats)
%
% SERVICES.RECONSTRUCTORMANAGER
%
%	Mock up of the RECONSTRUCTOR manager service
%
% INPUT
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'services.reconstructorManager';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

tTotal = tic();
fprintf(1, '\n STARTING RECONSTRUCTOR MANAGER on %s', nats.experimentName);
fprintf(1, '\n');

%% ALL INSTANCE SCALING AND CREATION
% here just execute the job added to the queue
[ok] = services.reconstructorInstance(nats);

%% REPORT OK
nats.simulatorService = ok;
fprintf(1, '\n');
fprintf(1, '\n RECONSTRUCTOR MANAGER on %s DONE', nats.experimentName);
fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
fprintf(1, '\n');
