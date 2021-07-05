%
% SEQUENCE PARAMETER TESTING
%
%	Runs the tests for the sequence package.
%   Test are defined in +sequence/+test folder
%   It will generate a code coverage report in ./testCoverage/sequence/
%
%   These test generate sequences for a bunch of parameters, and make 
%   sure that there is no crashing or error due to bad combination 
%   of input parameters.
%
%========================  CORSMED AB © 2020 ==============================

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
import matlab.unittest.plugins.codecoverage.CoverageReport

% create folders for report coverage
reportFolder = 'testCoverage';
if ~exist(reportFolder,'dir')
    mkdir(reportFolder);
end
if ~exist(sprintf('%s/sequence',reportFolder),'dir')
    mkdir(reportFolder,'sequence');
end
% define suite from a package
suite = TestSuite.fromPackage('sequence.test');
runner = TestRunner.withTextOutput;
% coverage plugin
genCobertura = 0; % generates report in cobertura format xml
if genCobertura
    reportFile      = './testCoverage/sequence/sequenceCoverageResults.xml';
    reportFormat    = CoberturaFormat(reportFile);
else
    reportFile      = 'sequenceCoverageResults.html';
    reportFormat    = CoverageReport('testCoverage/sequence', 'MainFile', reportFile );
end
reportPlugin    = CodeCoveragePlugin.forFolder('./+sequence', ...
    'IncludingSubfolders',true, 'Producing',reportFormat);
runner.addPlugin( reportPlugin );
% run
tTest = tic();
result = runner.run(suite);
tTest = toc(tTest);

fprintf(1, '\n Sequence Parameter Testing done');
fprintf(1, '\n    Elapsed Time  %g', tTest);
fprintf(1, '\n    Passed        %d', nnz([result.Passed]) );
fprintf(1, '\n    Failed        %d', nnz([result.Failed]) );
fprintf(1, '\n    Total Tests   %d', numel(result) );
fprintf(1, '\n\n');
