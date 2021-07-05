%
% EDUTOOL.REGRESSION
%
%	Runs the regression tests for the EduTool App
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
regressionFolder = '/efs-mount-point/S20/REGRESSION/eduTool';

%% Define the strings to parse:
%  the experiments containing any of the strings separated by + will be run
testString = 'defaultParamsVarys';

%% run the batch
tRegression = tic();
[regResult] = eduTool.test.runRegression( regressionFolder,testString, 'test');
tRegression = toc(tRegression);

%% report results
clc;
fprintf(1, '\n');
fprintf(1, '\n -------------------------------------------------------------------------------------------------');
fprintf(1, '\n  Regression results');
for ii = 1:regResult.numExperiments
    currentName = regResult.experiment{ii}.name;
    expName = sprintf('%75s', ' '); % for nice formatting
    expName(1:length(currentName)) = currentName(:);
    fprintf(1, '\n   %s  :  %s \t-- %s ',...
        regResult.experiment{ii}.result,...
        expName,...
        regResult.experiment{ii}.note);
    
    % print failing details if needed
    if regResult.experiment{ii}.report.pass == 0
        reportFields = fieldnames(regResult.experiment{ii}.report);
        % loop on the field names and call recursive for structs, otherwise assign
        for ff = 1:length(reportFields)
            if ~strcmpi(reportFields{ff},'pass')
                fprintf(1, '\n %30s  :  %s ',...
                    reportFields{ff}, regResult.experiment{ii}.report.(reportFields{ff}));
            end
        end
    end

end
fprintf(1, '\n -------------------------------------------------------------------------------------------------');
fprintf(1, '\n');
fprintf(1, '\n Regression run completed');
fprintf(1, '\n  Folder     : %s', regressionFolder);
fprintf(1, '\n  Query      : %s', testString);
fprintf(1, '\n  Passed     : %d/%d', regResult.numPass,regResult.numExperiments);
fprintf(1, '\n  Failed     : %d/%d', regResult.numFail,regResult.numExperiments);
fprintf(1, '\n  Total time : %.3fs', tRegression);
fprintf(1, '\n -------------------------------------------------------------------------------------------------');
fprintf(1, '\n');
