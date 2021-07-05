function [expInfo] = experimentExecutionV2(...
    anatomicalModel, coilSystem, mrSystem, ...
    acquisition, expControl, tagsStruct, parpoolConn)
%
% EDUTOOL.RUN.EXPERIMENTEXECUTION
%
%	Runs a experiment.
%
% INPUT
%   sessionData         solution struct with initial data
%   anatomicalModel     struct with the anatomical model
%   coilModel           struct with models of the different coils
%
% OUTPUT
%   expInfo             experiment info for the user
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'eduTool.run.experimentExecution';
if (nargin < 5)
    ME = MException('eduTool:wrongArgCount',...
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

%% SLICER
[spinModel, pulseSequence, motionModel, sarReport, imageData, reconData, ...
    expControl] = eduTool.run.slicer( anatomicalModel, coilSystem, mrSystem, ...
    acquisition, expControl );

%% check status
eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)

%% SIMULATOR
[simSignal] = eduTool.run.engine( spinModel, pulseSequence, motionModel,...
    expControl, tagsStruct, parpoolConn );

%% check status
eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)

%% RECONSTRUCTION
[reconData] = eduTool.run.recon(...
    simSignal, reconData, expControl );

%% check status
eduTool.frontend.checkExperimentStatus(expControl.connLocalDB,expControl.experimentID)

%% DICON GENERATION
[imageData] = eduTool.run.dicom(...
    reconData, imageData, expControl );

%% update progress bar to full
expControl.progress = 100;
eduTool.frontend.progressUpdate(expControl);

expInfo = 1;

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for experiment %d, elapsed time %.3fs',...
        functionName, expControl.experimentID, toc(tTotal));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
