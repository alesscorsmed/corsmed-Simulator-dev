function timeSolution = sbrEngine(spinModel,pulseSequence,motionModel,...
    expControl,GPUindex)
%
% SBR.RUN.ENGINE
%
%	Runs the simulator.
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sbr.run.engine';
if (nargin < 4)
    ME = MException('sbr:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% if there is no DB connection
if ~isfield(expControl,'connLocalDB')
    expControl.connLocalDB = [];
end

%% run the simulation
timeSolution    = [];
timeSolution    = simulator.runSimulation(spinModel, pulseSequence, ...
    motionModel, expControl, timeSolution, 1 );

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Number of Jobs    %d', simJobs.numJobs);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end