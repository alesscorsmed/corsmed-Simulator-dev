function [results] = runUnitTest(type)
%
% SPINTWIN..RUNUNITTEST
%
%
%   
%
%========================  CORSMED AB © 2020 ==============================
%

functionName = 'spinTwin.runUnitTest';
if (nargin < 1 || isempty(type))
    type='all';
end

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport


%% prepare empty suites
fidSuite                = [];
flipAngleSuite          = [];
sliceSelectionSuite     = [];
signalIntegrationSuite  = [];
phaseRefocusSuite       = [];
diffusionSuite          = [];

%% add all bloch related testing
if strcmpi(type,'all') || strcmpi(type,'fid')
    fidSuite = TestSuite.fromFile(...
        './+spinTwin/+test/+unit/fidTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'flipAngle')
    flipAngleSuite = TestSuite.fromFile(...
        './+spinTwin/+test/+unit/flipAngleTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'sliceSelection')
    sliceSelectionSuite = TestSuite.fromFile(...
        './+spinTwin/+test/+unit/sliceSelectionTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'phaseRefocus')
    phaseRefocusSuite = TestSuite.fromFile(...
        './+spinTwin/+test/+unit/phaseRefocusTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'signalIntegration')
    signalIntegrationSuite = TestSuite.fromFile(...
        './+spinTwin/+test/+unit/signalIntegrationTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'diffusion')
    diffusionSuite = TestSuite.fromFile(...
        './+spinTwin/+test/+unit/pgseDiffusionTest.m');
end

%% join suits, run the tests and present table
simulatorSuite = [ fidSuite, flipAngleSuite, sliceSelectionSuite, ...
    signalIntegrationSuite, phaseRefocusSuite, diffusionSuite ];

%% prepare the runner and add the coverage
runner = TestRunner.withNoPlugins;
runner.addPlugin(CodeCoveragePlugin.forFolder('+spinTwin', ...
    'IncludingSubfolders',true, 'Producing', ...
    CoverageReport('+spinTwin/testCoverage', 'MainFile','spinTwin.html')));

%% run the suite
results = runner.run(simulatorSuite);

% final 
if any([results.Failed] == true)
    table(results)
    fprintf(1, '\n OH NO! %d fails out of %d tests... check the results above\n\n', ...
        nnz([results.Failed]), length(results));
else
    fprintf(1, '\n AWESOME! all %d tests passed... good job! \n\n', ...
        nnz([results.Passed]));
end

