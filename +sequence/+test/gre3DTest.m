%
%	Testing to verify generation of sequence under a variation of
%	parameters.
%
%   Uses parameterized class test
%   for more info: 
%   https://www.mathworks.com/help/matlab/matlab_prog/create-basic-parameterized-test.html
%
%   
%
%========================  CORSMED AB © 2020 ==============================
%

% define the class
classdef gre3DTest < matlab.unittest.TestCase
    
    % parameterization of the tests
    properties (TestParameter)
        %% we will run sequential mode, so each column will be a test
        fovFE               = {0.15,        0.15,               0.20,               0.20,       0.23,       0.30,               0.30}; % field of view in m
        fovPE               = {0.15,        0.20,               0.20,               0.20,       0.20,       0.25,               0.30}; % field of view in m
        numFE               = { 128,        256,                128,                256,        256,        129,                256}; % freq encodings
        numPE               = { 128,        128,                128,                256,        255,        196,                128};  % phase encodings
        sliceThickness      = {6e-3,        6e-3,               5e-3,               3e-3,       5e-3,       9e-3,               6e-3};
        sliceGap            = {2e-3,        1e-3,               1e-3,               3e-3,       1e-3,       1e-3,               3e-3};
        numSlices           = {   1,        5,                  50,                 10,         25,         5,                  15};
        rxBW                = {50e3,        75e3,               100e3,              50e3,       250e3,      150e3,              200e3 };
        TR                  = {6.5e-3,      15e-3,              7e-3,               25e-3,      50e-3,      150e-3,             400e-3};
        TE                  = {3.5e-3,      10e-3,              4e-3,               15e-3,      25e-3,      5e-3,               7e-3};
        foldoverSuppr       = {'no',        'yes',              'no',               'no',       'yes',      'yes',              'no'};
        parallelImaging     = {'no',        'yes',              'yes',              'yes',      'no',       'no',               'no'};
        rFactor             = {1,           1,                  2,                  4,          1,          3,                  2 };
        partialFourier      = {'none',      'phaseConjugate',   'readConjugate',    'none',     'none',     'phaseConjugate',   'readConjugate'};
        fFactor             = {1.0,         0.8,                0.6,                1.0,        1.0,        0.8,                0.6};
    end
    
    % Functions with unit tests
    methods (Test, ParameterCombination = 'sequential')
        
        %% test GRE generation with Edge cases
        function runGre3DEdgeTest(testCase,...
                fovFE,...
                fovPE,...
                numFE,...
                numPE,...
            	sliceThickness,...
                sliceGap,...
            	numSlices,...
                rxBW,...
                TR,...
                TE,...
                foldoverSuppr,...
                parallelImaging,...
                rFactor,...
                partialFourier,...
                fFactor)
            
            %% initialize
            expControl          = data.experiment.initializeExpControl();
            acquisition         = data.experiment.initializeAcquisition();
            % fixed vars
            acquisition.data.pulseSeqFamilyName = 'GRE-3D';
            acquisition.data.pulseSeqActualName = 'GRE-3D';
            acquisition.data.is3D               = 1;
            % 
            expControl.debug.devMode = 0;
            expControl.debug.debugMode = 0;
            
            %% general info for the test
            acquisition.data.fovFE              = fovFE;
            acquisition.data.fovPE              = fovPE;
            acquisition.data.numFE              = numFE;
            acquisition.data.numPE              = numPE;
            
            acquisition.data.sliceGap           = sliceGap;
            acquisition.data.sliceThickness     = sliceThickness;
            acquisition.data.numSlices          = numSlices;
            
            acquisition.data.fovSE = acquisition.data.numSlices ...
                * acquisition.data.sliceThickness;
            acquisition.data.numSE = acquisition.data.numSlices;
            
            acquisition.data.TR                 = TR;
            acquisition.data.TE                 = TE;
            acquisition.data.rxBW               = rxBW;
            
            acquisition.data.foldoverSuppr      = foldoverSuppr;
            acquisition.data.parallelImaging    = parallelImaging;
            acquisition.data.rFactor            = rFactor;
            acquisition.data.partialFourier     = partialFourier;
            acquisition.data.fFactor            = fFactor;

            try
                %% initialize the encoding
                [encodingData] = encoder.generateEncodingPlan( ...
                    acquisition, expControl );
                %% generate the sequence
                [~] = sequence.generatePulseSequence( ...
                    acquisition, encodingData, expControl.mrSystem, expControl );
                %%
                pass = true;
            catch ME
                fprintf(1, '\n ERROR: %s (%d) : %s', ME.stack(1).file, ME.stack(1).line, ME.message);
                pass = false;
            end

            %% verify passing test
            testCase.verifyTrue(pass);
            
        end

        
        %% test GRE generation with General cases
        function runGre3DGeneralTest(testCase)
            
            %% initialize
            expControl          = data.experiment.initializeExpControl();
            acquisition         = data.experiment.initializeAcquisition();
            % fixed vars
            acquisition.data.pulseSeqFamilyName = 'GRE-3D';
            acquisition.data.pulseSeqActualName = 'GRE-3D';
            acquisition.data.is3D               = 1;
            % 
            expControl.debug.devMode = 0;
            expControl.debug.debugMode = 0;
            
            %% general info for the test            
            acquisition.data.sliceGap           = 3e-3;
            
            acquisition.data.foldoverSuppr      = 'no';
            acquisition.data.parallelImaging    = 'no';
            acquisition.data.rFactor            = 1;
            acquisition.data.partialFourier     = 'none';
            acquisition.data.fFactor            = 1.0;
            
            numPass = 0;
            numFail = 0;
            numWarn = 0;
            fprintf(1,'\n');
            
            %% loop on general cases
            for lfovFE = linspace(0.11, 0.42, 3)
            for lfovPE = linspace(0.15, 0.35, 3)
            for lnumFE = [96 128 256]
            for lnumPE = [64 128]
            for lTE = linspace(4e-3, 30e-3, 4)
            for lTR = linspace(35e-3, 4, 4)
            for lrxBW = [75 200]*1e3
            for lsliceThickness = [3 6]*1e-3
            for lnumSlices = [3 6 10]
            
            %% general info for the test
            acquisition.data.fovFE              = lfovFE;
            acquisition.data.fovPE              = lfovPE;
            acquisition.data.numFE              = lnumFE;
            acquisition.data.numPE              = lnumPE;
            acquisition.data.TR                 = lTR;
            acquisition.data.TE                 = lTE;
            acquisition.data.rxBW               = lrxBW;
            acquisition.data.numSlices          = lnumSlices;
            acquisition.data.sliceThickness     = lsliceThickness;
            acquisition.data.fovSE              = lnumSlices*lsliceThickness;
            acquisition.data.numSE              = lnumSlices;

            try
                %% initialize the encoding
                [encodingData] = encoder.generateEncodingPlan( ...
                    acquisition, expControl );
                %% generate the sequence
                [~] = sequence.generatePulseSequence( ...
                    acquisition, encodingData, expControl.mrSystem, expControl );
                %%
                numPass = numPass + 1;
            catch ME
                %% check if error is user warning
                errorInfo= strsplit(ME.identifier,':');
                if strcmpi(errorInfo{1}, 'userwarning')
                    numWarn = numWarn + 1;
                else
                    %% otherwise error
                    fprintf(1, '\n GRE ERROR: %s (%d) : %s', ME.stack(1).file, ME.stack(1).line, ME.message);
                    fprintf(1, ['\n    PARAMETERS:',...
                        '  fovFE=%f  fovPE=%f  numFE=%d   numPE=%d',...
                        '  TE=%.3g  TR=%.3g  rxBW=%.2g  sliceThicness=%.1g slices=%d'],...
                        lfovFE, lfovPE, lnumFE, lnumPE, lTE, lTR, lrxBW, lsliceThickness, lnumSlices);
                    numFail = numFail + 1;
                end
            end
            
            numTotal = numPass+numFail+numWarn;
            if rem(numTotal,100) == 0
                fprintf(1,'.');
            end
            if rem(numTotal,1000) == 0
                fprintf(1,'  ');
            end
            if rem(numTotal,10000) == 0
                fprintf(1,'\n');
            end
            
            end
            end
            end
            end
            end
            end
            end
            end
            end

            fprintf(1, '\n TOTAL %d : PASS %d / WARN %d / FAIL %d \n',...
                numPass+numFail+numWarn, numPass, numWarn, numFail);
            %% verify number of Fails = 0
            testCase.verifyEqual(numFail, 0);
            
        end

        
    end
end