%
% SPIN TWIN TESTING
%
%	Runs the tests for the SpinTwin package.
%   Test are defined in +spinTwin/+test/+unit folder
%   It will generate a code coverage report in ./testCoverage/spinTwin/
%
%   These are physical-based test that verify the CUDA kernels.
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
if ~exist(sprintf('%s/spinTwin',reportFolder),'dir')
    mkdir(reportFolder,'spinTwin');
end
% define suite from a package
suite = TestSuite.fromPackage('spinTwin.test.unit');
runner = TestRunner.withTextOutput;
% coverage plugin
genCobertura = 0; % generates report in cobertura format xml
if genCobertura
    reportFile      = './testCoverage/spinTwin/spinTwinCoverageResults.xml';
    reportFormat    = CoberturaFormat(reportFile);
else
    reportFile      = 'spinTwinCoverageResults.html';
    reportFormat    = CoverageReport('testCoverage/spinTwin', 'MainFile', reportFile );
end
reportPlugin    = CodeCoveragePlugin.forFolder('./+spinTwin', ...
    'IncludingSubfolders',true, 'Producing',reportFormat);
runner.addPlugin( reportPlugin );
% run
result = runner.run(suite);