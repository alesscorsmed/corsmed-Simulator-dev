function [results] = runTestSuite(type)
%
% SIMULATOR.TEST.RUNTESTSUITE
%
%
%   
%
%========================  CORSMED AB © 2020 ==============================
%

functionName = 'simulator.test.runTestSuite';
if (nargin < 1 || isempty(type))
    type='all';
end

%% prepare empty suites
basicBlochSuite         = [];
flipAngleSuite          = [];
sliceSelectionSuite     = [];
signalIntegrationSuite  = [];
refocusingSuite         = [];
diffusionSuite          = [];
crossDiffusionSuite     = [];

%% add all bloch related testing
if strcmpi(type,'all') || strcmpi(type,'flipAngle')
    flipAngleSuite = matlab.unittest.TestSuite.fromFile(...
        './+simulator/+bloch/+test/blochFlipAngleTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'sliceSelection')
    sliceSelectionSuite = matlab.unittest.TestSuite.fromFile(...
        './+simulator/+bloch/+test/blochSliceSelectionTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'refocusing')
    refocusingSuite = matlab.unittest.TestSuite.fromFile(...
        './+simulator/+bloch/+test/blochRefocusingTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'signalIntegration')
    signalIntegrationSuite = matlab.unittest.TestSuite.fromFile(...
        './+simulator/+bloch/+test/signalIntegrationTest.m');
end


%% add all diffusion related testing
if strcmpi(type,'all') || strcmpi(type,'diffusion')
    crossDiffusionSuite = matlab.unittest.TestSuite.fromFile(...
        './+simulator/+diffusion/+test/crossComparisonTest.m');
end
if strcmpi(type,'all') || strcmpi(type,'diffusion')
    diffusionSuite = matlab.unittest.TestSuite.fromFile(...
        './+simulator/+diffusion/+test/pgseDiffusionTest.m');
end

%% join suits, run the tests and present table
simulatorSuite = [ basicBlochSuite, flipAngleSuite, sliceSelectionSuite, ...
    signalIntegrationSuite, refocusingSuite, crossDiffusionSuite, diffusionSuite ];
results = simulatorSuite.run;
% final 
if any([results.Failed] == true)
    table(results)
    fprintf(1, '\n OH NO! %d fails out of %d tests... check the results above\n\n', ...
        nnz([results.Failed]), length(results));
else
    fprintf(1, '\n AWESOME! all %d tests passed... good job! \n\n', ...
        nnz([results.Passed]));
end

