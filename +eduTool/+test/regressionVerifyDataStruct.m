function [testControl] = regressionVerifyDataStruct( targetStruct, testControl, structName )
%
% EDUTOOL.TEST.REGRESSIONVERIFYSTRUCT
%
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%

%% Load reference data from file
fileName = sprintf('%s/golden/%s_%s.mat', ...
    testControl.regFolder, testControl.testName, structName);
try
    load(fileName, 'dataStruct');
    doTest = 1;
catch
    testControl.status.(structName) = sprintf('FAIL : unable to load reference %s', fileName);
    testControl.status.pass = 0;
    doTest = 0;
end

if doTest
    %% compare both structs
    throwError  = 1; % if there are differences, throw
    relTol      = 1e-4; % relative error for the norm
    try
        tools.misc.deepCompareDataStruct( ...
            targetStruct, dataStruct.(structName), relTol, throwError );
        % identical structs: pass regression
        testControl.status.(structName) = 'PASS';
    catch ME
        % fail regression
        testControl.status.(structName) = sprintf('FAIL : %s', ME.message);
        testControl.status.pass = 0;
    end
end
