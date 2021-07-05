%
% EDUTOOL.BATCHMAIN
%
%	Runs the Batch procedure for the EduTool App
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

%% Define the folder of the experiments
batchFolder = '/efs-mount-point/S20/TESTS/edutool';

%% Define the strings to parse:
%  the experiments containing any of the strings separated by + will be run
batchString = 'GRE+SAG';

%% Define if you want to update the json files
%  useful in case where new fields are added to expControl or acquisition
updateJson = 0;

%% run the batch
tBatch = tic();
[batchResult] = eduTool.batch.runBatch(batchFolder,batchString,updateJson);
tBatch = toc(tBatch);

%% report results
fprintf(1, '\n');
fprintf(1, '\n');
fprintf(1, '\n -----------------------------------------------------------------');
fprintf(1, '\n Batch run completed');
fprintf(1, '\n  Folder     : %s', batchFolder);
fprintf(1, '\n  Query      : %s', batchString);
fprintf(1, '\n  Passed     : %d/%d', batchResult.numPass,batchResult.numExperiments);
fprintf(1, '\n  Failed     : %d/%d', batchResult.numFail,batchResult.numExperiments);
fprintf(1, '\n  Total time : %.3fs', tBatch);
fprintf(1, '\n -----------------------------------------------------------------');
fprintf(1, '\n  Results');
for ii = 1:batchResult.numExperiments
    currentName = batchResult.experiment{ii}.name;
    expName = sprintf('%50s', ' '); % for nice formatting
    expName(1:length(currentName)) = currentName(:);
    fprintf(1, '\n   %s  :  %s \t-- %s ',...
        batchResult.experiment{ii}.result,...
        expName,...
        batchResult.experiment{ii}.note);
end
fprintf(1, '\n -----------------------------------------------------------------');
fprintf(1, '\n');


