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
classdef ssfseTest < matlab.unittest.TestCase
    
    % parameterization of the tests
    properties (TestParameter)
        %% we will run sequential mode, so each column will be a test
        fovFE               = {0.15,        0.15,               0.20,               0.20,       0.23,       0.30,               0.30}; % field of view in m
        fovPE               = {0.15,        0.20,               0.20,               0.20,       0.20,       0.25,               0.30}; % field of view in m
        numFE               = { 128,        256,                128,                256,        256,        129,                256}; % freq encodings
        numPE               = { 128,        129,                129,                256,        255,        196,                128};  % phase encodings
        sliceThickness      = {6e-3,        6e-3,               5e-3,               3e-3,       5e-3,       9e-3,               6e-3};
        sliceGap            = {2e-3,        1e-3,               1e-3,               3e-3,       1e-3,       1e-3,               3e-3};
        numSlices           = {   1,        5,                  50,                 10,         25,         5,                  15};
        rxBW                = {50e3,        75e3,               100e3,              50e3,       250e3,      150e3,              200e3 };
        TR                  = {900e-3,      1200e-3,            4000e-3,            3500e-3,    3000e-3,    4500e-3,            800e-3};
        TE                  = {7.0e-3,      10e-3,              7e-3,               10e-3,      8.0e-3,     25e-3,              10e-3};
        foldoverSuppr       = {'no',        'yes',              'no',               'no',       'yes',      'yes',              'no'};
        partialFourier      = {'none',      'phaseConjugate',   'phaseConjugate',   'none',     'none',     'phaseConjugate',   'phaseConjugate'};
        fFactor             = {1.0,         0.8,                0.6,                1.0,        1.0,        0.8,                0.55};
        fatsat              = {'none',      'none',             'none',             'lipidir',  'lipidir',  'none',             'none'};
        TI                  = {215e-3,      215e-3,             215e-3,             215e-3,      500e-3,    215e-3,             215e-3};

    end
    
    % Functions with unit tests
    methods (Test, ParameterCombination = 'sequential')
        
        %% test TSE generation with Edge cases
        function runSSFSEEdgeTest(testCase,...
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
                partialFourier,...
                fFactor,...
                fatsat,...
                TI)
            
            %% initialize
            expControl          = data.experiment.initializeExpControl();
            acquisition         = data.experiment.initializeAcquisition();
            % fixed vars
            acquisition.data.pulseSeqFamilyName = 'SS-FSE';
            acquisition.data.pulseSeqActualName = 'SS-FSE';
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
            acquisition.data.partialFourier     = partialFourier;
            acquisition.data.fFactor            = fFactor;
            
            acquisition.data.fatsat             = fatsat;
            if strcmpi(fatsat,'lipidir')
                acquisition.prepIR.apply        = 1;
                acquisition.prepIR.TI           = TI;
            end
            
            acquisition.data.shotTR             = TR;

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

        
        %% test SSFSE generation with General cases
        function runSSFSEGeneralTest(testCase)
            
            %% initialize
            expControl          = data.experiment.initializeExpControl();
            acquisition         = data.experiment.initializeAcquisition();
            % fixed vars
            acquisition.data.pulseSeqFamilyName = 'SS-FSE';
            acquisition.data.pulseSeqActualName = 'SS-FSE';
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
            
            numPass = 0;
            numFail = 0;
            numWarn = 0;
            fprintf(1,'\n');
            
            %% loop on general cases
            for lfovFE = linspace(0.15, 0.35, 3)
            for lfovPE = linspace(0.15, 0.35, 3)
            for lnumFE = [129 256 512]
            for lnumPE = [96  129]
            for lTE = linspace(7e-3, 25e-3, 3)
            for lTR = [1.5, 4.0]
            for lrxBW = [75 250]*1e3
            for lsliceThickness = 6e-3
            for lTI = [0, 215, 500]*1e-3
            for lpartialFourier = {'none', 'phaseConjugate'}
            for lfFactor = [0.6, 0.8, 1.0]

            
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

            acquisition.data.shotTR             = lTR;

            acquisition.data.partialFourier     = lpartialFourier{1};
            acquisition.data.fFactor            = lfFactor;
            

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
                    fprintf(1, '\n SS-FSE ERROR: %s (%d) : %s', ME.stack(1).file, ME.stack(1).line, ME.message);
                    fprintf(1, ['\n    PARAMETERS:',...
                        '  fovFE=%f  fovPE=%f  numFE=%d   numPE=%d',...
                        '  TE=%.3g  TR=%.3g  rxBW=%.2g  sliceThicness=%.1g',...
                        '  FATSAT=%s  TI=%g  PartialFourier=%s, FourierFactor'],...
                        lfovFE, lfovPE, lnumFE, lnumPE, lTE, lTR, lrxBW, lsliceThickness,...
                        acquisition.data.fatsat, lTI, lpartialFourier{1}, lfFactor);
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
            
            %% verify number of Fails = 0
            fprintf(1, '\n TOTAL %d : PASS %d / WARN %d / FAIL %d \n',...
                numPass+numFail+numWarn, numPass, numWarn, numFail);
            testCase.verifyEqual(numFail, 0);
            
        end

        
    end
end