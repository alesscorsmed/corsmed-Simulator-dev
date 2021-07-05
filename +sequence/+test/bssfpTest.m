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
classdef bssfpTest < matlab.unittest.TestCase
    
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
        TR                  = {4e-3,        15e-3,              2.72e-3,            25e-3,      50e-3,      150e-3,             400e-3};
        TE                  = {2e-3,        10e-3,              1.36e-3,            15e-3,      25e-3,      5e-3,               7e-3};
        foldoverSuppr       = {'no',        'yes',              'no',               'no',       'yes',      'yes',              'no'};
        parallelImaging     = {'no',        'yes',              'no',               'yes',      'no',       'no',               'no'};
        rFactor             = {1,           1,                  2,                  4,          1,          3,                  2 };
        partialFourier      = {'none',      'phaseConjugate',   'none',             'none',     'none',     'phaseConjugate',   'readConjugate'};
        fFactor             = {1.0,         0.8,                0.6,                1.0,        1.0,        0.8,                0.6};
        fatsat              = {'none',      'none',             'none',             'lipidir',  'lipidir',  'none',             'none'};
        TI                  = {215e-3,      215e-3,             215e-3,             215e-3,      500e-3,    215e-3,             215e-3};
        prepBSSFP           = {'ramp',      'half',             'ramp',             'half',      'ramp',    'half',             'ramp'};
        dummySSFP           = {3,           5,                  10,                  20,          10,         8,                 6 };

    end
    
    % Functions with unit tests
    methods (Test, ParameterCombination = 'sequential')
        
        %% test BSSFP generation with Edge cases
        function runBSSFPEdgeTest(testCase,...
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
                fFactor,...
                fatsat,...
                TI,...
                prepBSSFP,...
                dummySSFP)
            
            %% initialize
            expControl          = data.experiment.initializeExpControl();
            acquisition         = data.experiment.initializeAcquisition();
            % fixed vars
            acquisition.data.pulseSeqFamilyName = 'bSSFP';
            acquisition.data.pulseSeqActualName = 'bSSFP';
            acquisition.mainRF.duration         = 0.5e-3; % fixed for bssfp
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
            
            acquisition.data.fovSE = acquisition.data.sliceThickness;
            acquisition.data.numSE = 0;
            
            acquisition.data.TR                 = TR;
            acquisition.data.TE                 = TE;
            acquisition.data.rxBW               = rxBW;
            
            acquisition.data.foldoverSuppr      = foldoverSuppr;
            acquisition.data.parallelImaging    = parallelImaging;
            acquisition.data.rFactor            = rFactor;
            acquisition.data.partialFourier     = partialFourier;
            acquisition.data.fFactor            = fFactor;
            
            acquisition.data.fatsat             = fatsat;
            if strcmpi(fatsat,'lipidir')
                acquisition.prepIR.apply        = 1;
                acquisition.prepIR.TI           = TI;
            end
            acquisition.data.prepBSSFP          = prepBSSFP;
            acquisition.data.dummySSFP          = dummySSFP;

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
        function runBSSFPGeneralTest(testCase)
            
            %% initialize
            expControl          = data.experiment.initializeExpControl();
            acquisition         = data.experiment.initializeAcquisition();
            % fixed vars
            acquisition.data.pulseSeqFamilyName = 'bSSFP';
            acquisition.data.pulseSeqActualName = 'bSSFP';
            acquisition.mainRF.duration         = 0.5e-3; % fixed for bssfp
            % 
            expControl.debug.devMode = 0;
            expControl.debug.debugMode = 0;
            
            %% general info for the test            
            acquisition.data.sliceGap           = 3e-3;
            acquisition.data.numSlices          = 1;
            acquisition.data.numSE = 0;    
            
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
            for lfovFE = linspace(0.15, 0.4, 3)
            for lfovPE = linspace(0.15, 0.4, 3)
            for lnumFE = [128 256 512]
            for lnumPE = [128 192]
            for lTE = linspace(1.3e-3, 12e-3, 3)
            for lTR = linspace(2.6e-3, 50e-3, 3)
            for lrxBW = [50 100 250]*1e3
            for lsliceThickness = [3 6]*1e-3
            for lTI = [0, 215, 800]*1e-3
            for lprepBSSFP = {'ramp', 'half'}
            for ldummySSFP = [5 10]
            
            %% general info for the test
            acquisition.data.fovFE              = lfovFE;
            acquisition.data.fovPE              = lfovPE;
            acquisition.data.numFE              = lnumFE;
            acquisition.data.numPE              = lnumPE;
            acquisition.data.TR                 = lTR;
            acquisition.data.TE                 = lTE;
            acquisition.data.rxBW               = lrxBW;
            acquisition.data.sliceThickness     = lsliceThickness;
            acquisition.data.fovSE = acquisition.data.sliceThickness;
            
            if lTI > 0
                acquisition.data.fatsat         = 'lipidir';
                acquisition.prepIR.apply        = 1;
            else
                acquisition.data.fatsat         = 'none';
            end
            acquisition.data.prepBSSFP          = lprepBSSFP{1};
            acquisition.data.dummySSFP          = ldummySSFP;

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
                        '  TE=%.3g  TR=%.3g  rxBW=%.2g  sliceThicness=%.1g',...
                        '  FATSAT=%s  TI=%g  prepBSSFP=%s  dummySSFP=%d'],...
                        lfovFE, lfovPE, lnumFE, lnumPE, lTE, lTR, lrxBW, lsliceThickness,...
                        acquisition.data.fatsat, lTI, lprepBSSFP{1}, ldummySSFP);
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
            end
            end
            
            fprintf(1, '\n TOTAL %d : PASS %d / WARN %d / FAIL %d \n',...
                numPass+numFail+numWarn, numPass, numWarn, numFail);
            %% verify number of Fails = 0
            testCase.verifyEqual(numFail, 0);
            
        end

        
    end
end