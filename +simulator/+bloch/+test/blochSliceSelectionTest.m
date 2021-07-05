%
% SIMULATOR.BLOCH.TEST.BLOCHSLICESELECTIONTEST
%
%	Functional unit testing for the Bloch Simulator.
%   Verifies correctness of RF flipping with slice selection
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
classdef blochSliceSelectionTest < matlab.unittest.TestCase   
    
    % parameterization of the tests
    properties (TestParameter)
        expControl      = {data.expControl.initialize()};
        t1Value         = {10.0};
        t2Value         = {7.0};
        gradAmp         = {0.03}
        gradSlew        = {0, 150};
        sliceThickness  = {4e-3, 6e-3};
        rfDuration      = {3e-3, 5e-3};
        rfFlipAngle     = {90};
        rfPhase         = {0.0, pi/2};
        rfCycles        = {2, 6, 10};       
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test the the flip angle without slice selection
        %% ANALYTICAL kernel
        function runAnaliticalSliceSelectionTest(testCase,expControl,...
                t1Value,t2Value,gradAmp,gradSlew,...
                rfDuration,rfFlipAngle,rfPhase,...
                rfCycles, sliceThickness)
            
            expControl.debug.debugMode    = 0;
            motionModel = [];
            
            %% mr System specs
            [mrSystem] = data.mrSystem.initialize();
            mrSystem.maxGStrenght   = gradAmp;
            mrSystem.SlewRate       = gradSlew;
            
            %% generate the sequence
            mainRF      = data.acquisition.initializeRF(rfFlipAngle);
            %% assign values
            mainRF.flipAngle          = rfFlipAngle; % in degrees
            mainRF.phase              = rfPhase; % in rads
            mainRF.duration           = rfDuration;
            mainRF.cycles             = rfCycles;
            mainRF.sliceThickness     = sliceThickness;
            mainRF.doSliceSelect      = 1;
            mainRF.preRWScale         = 0;
            mainRF.postRWScale        = 0;
            %% populate waveforms
            [pulseSequence]=sequence.dummy.singleRF(...
                mainRF, mrSystem, expControl);
            
            %% initial magnetization state
            b0                  = 1.0;
            mu                  = 1;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 0.0;
            
            %% domain
            fovX = 0.001;
            fovY = 0.001;
            fovZ = 2*sliceThickness;
            dx   = 1e-3;
            dy   = 1e-3;
            dz   = 1e-4;
            
            %% generate spin Model
            [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ, dx, dy, dz,...
                t1Value, t2Value, b0, mu, expControl.debug.debugMode);
            
            %% initial magnetizations
            % get sizes
            numCoils = spinModel.slice{1}.model.numCoils;
            numIso   = spinModel.slice{1}.model.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
                       
            %% call the standard analytical kernel
            solAnalytical = data.solution.initialize(numIso,numRxs,numCoils);
            solAnalytical.indexes    = [];
            solAnalytical.Mm         = 1*M0m;
            solAnalytical.Mp         = 1*M0p;
            solAnalytical.Mx         = 1*M0x;
            solAnalytical.My         = 1*M0y;
            solAnalytical.Mz         = 1*M0z;
            solAnalytical.dMpDx      = 0*M0p;
            solAnalytical.dMpDy      = 0*M0p;
            solAnalytical.dMpDz      = 0*M0p;
            solAnalytical.Sx         = zeros(numCoils*numRxs,1);
            solAnalytical.Sy         = zeros(numCoils*numRxs,1);
            solAnalytical.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solAnalytical] = simulator.bloch.kernel.runAnalytical(...
                solAnalytical, pulseSequence, spinModel.slice{1}.model,...
                motionModel, expControl);

            %% target magnetizations
            refMm = initialMagnitude*sin(rfFlipAngle*pi/180)*ones(numIso,1);
            
            %% depending on the cycles, better or worst selectivity
            switch rfCycles
                case 10
                    % within 5% of target FA in 95% of the slice thickness
                    allowedError        = 0.05; 
                    sliceSelectivity    = 0.95; 
                    % magnitude reduction off slice is > 50%
                    offSliceReduction   = 0.50; 
                case 6
                    % within 5% of target FA in 90% of the slice thickness
                    allowedError        = 0.05; 
                    sliceSelectivity    = 0.90; 
                    % magnitude reduction off slice is > 50%
                    offSliceReduction   = 0.50; 
                    
                otherwise
                    % within 10% of target FA in 80% of the slice thickness
                    allowedError        = 0.10; 
                    sliceSelectivity    = 0.80; 
                    % magnitude reduction off slice is > 50%
                    offSliceReduction   = 0.50;
            end
            
            %% verify right slice selection
            idxInSlice = find(abs(spinModel.slice{1}.model.z) < ...
                sliceSelectivity*sliceThickness/2);
            errorInSlice = abs(solAnalytical.Mm(idxInSlice)-refMm(idxInSlice));
            testCase.verifyLessThan(errorInSlice,...
                allowedError*max(abs(refMm)));
            %% verify right off slice dampening
            idxOffSlice = find(abs(spinModel.slice{1}.model.z) > ...
                sliceThickness/2);
            testCase.verifyLessThan(solAnalytical.Mm(idxOffSlice),...
                offSliceReduction*max(abs(refMm)));
            
        end
        
        %% test the the flip angle without slice selection
        %% PHASOR kernel
        function runPhasorSliceSelectionTest(testCase,expControl,...
                t1Value,t2Value,gradAmp,gradSlew,...
                rfDuration,rfFlipAngle,rfPhase,rfCycles,sliceThickness)
            
            expControl.debug.debugMode    = 0;
            motionModel = [];
            
            %% mr System specs
            [mrSystem] = data.mrSystem.initialize();
            mrSystem.maxGStrenght   = gradAmp;
            mrSystem.SlewRate       = gradSlew;
            
            %% generate the sequence
            mainRF      = data.acquisition.initializeRF(rfFlipAngle);
            %% assign values
            mainRF.flipAngle          = rfFlipAngle; % in degrees
            mainRF.phase              = rfPhase; % in rads
            mainRF.duration           = rfDuration;
            mainRF.cycles             = rfCycles;
            mainRF.sliceThickness     = sliceThickness;
            mainRF.doSliceSelect      = 1;
            mainRF.preRWScale         = 0;
            mainRF.postRWScale        = 0;
            %% populate waveforms
            [pulseSequence]=sequence.dummy.singleRF(...
                mainRF, mrSystem, expControl);
            
            %% initial magnetization state
            b0                  = 1.0;
            mu                  = 1;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 0.0;
            
            %% domain
            fovX = 0.001;
            fovY = 0.001;
            fovZ = 2*sliceThickness;
            dx   = 1e-3;
            dy   = 1e-3;
            dz   = 1e-4;
            
            %% generate spin Model
            [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ, dx, dy, dz,...
                t1Value, t2Value, b0, mu, expControl.debug.debugMode);
            
            %% initial magnetizations
            % get sizes
            numCoils = spinModel.slice{1}.model.numCoils;
            numIso   = spinModel.slice{1}.model.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
                       
            %% call the standard analytical kernel
            solPhasor = data.solution.initialize(numIso,numRxs,numCoils);
            solPhasor.indexes    = [];
            solPhasor.Mm         = 1*M0m;
            solPhasor.Mp         = 1*M0p;
            solPhasor.Mx         = 1*M0x;
            solPhasor.My         = 1*M0y;
            solPhasor.Mz         = 1*M0z;
            solPhasor.dMpDx      = 0*M0p;
            solPhasor.dMpDy      = 0*M0p;
            solPhasor.dMpDz      = 0*M0p;
            solPhasor.Sx         = zeros(numCoils*numRxs,1);
            solPhasor.Sy         = zeros(numCoils*numRxs,1);
            solPhasor.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solPhasor] = simulator.bloch.kernel.runPhasor(...
                solPhasor, pulseSequence,spinModel.slice{1}.model, ...
                motionModel, expControl);

            %% target magnetizations
            refMm = initialMagnitude*sin(rfFlipAngle*pi/180)*ones(numIso,1);
            
            %% depending on the cycles, better or worst selectivity
            switch rfCycles
                case 10
                    % within 5% of target FA in 95% of the slice thickness
                    allowedError        = 0.05; 
                    sliceSelectivity    = 0.95; 
                    % magnitude reduction off slice is > 50%
                    offSliceReduction   = 0.50; 
                case 6
                    % within 5% of target FA in 90% of the slice thickness
                    allowedError        = 0.05; 
                    sliceSelectivity    = 0.90; 
                    % magnitude reduction off slice is > 50%
                    offSliceReduction   = 0.50; 
                    
                otherwise
                    % within 10% of target FA in 80% of the slice thickness
                    allowedError        = 0.10; 
                    sliceSelectivity    = 0.80; 
                    % magnitude reduction off slice is > 50%
                    offSliceReduction   = 0.50;
            end
            
            %% verify right slice selection
            idxInSlice = find(abs(spinModel.slice{1}.model.z) < ...
                sliceSelectivity*sliceThickness/2);
            errorInSlice = abs(solPhasor.Mm(idxInSlice)-refMm(idxInSlice));
            testCase.verifyLessThan(errorInSlice,...
                allowedError*max(abs(refMm)));
            %% verify right off slice dampening
            idxOffSlice = find(abs(spinModel.slice{1}.model.z) > ...
                sliceThickness/2);
            testCase.verifyLessThan(solPhasor.Mm(idxOffSlice),...
                offSliceReduction*max(abs(refMm)));
            
        end
        
    end
end