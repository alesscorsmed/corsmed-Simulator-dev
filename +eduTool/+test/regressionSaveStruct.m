function [testControl] = regressionSaveStruct( targetStruct, testControl, structName )
%
% EDUTOOL.TEST.REGRESSIONSAVESTRUCT
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

% add version used in generation, and time stamp
dataStruct.version      = testControl.version;
dataStruct.timeStamp    = testControl.timeStamp;
dataStruct.(structName) = targetStruct;
try
    save(fileName, 'dataStruct', '-v7.3');
    testControl.status.(structName) = sprintf('PASS : saved reference %s', fileName);
catch
    testControl.status.(structName) = sprintf('FAIL : unable to save reference %s', fileName);
    testControl.status.pass = 0;
end