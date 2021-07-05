%
% SIMULATOR.DIFUSSION.TEST.CROSSCOMPARISONEST
%
%	Compares the results of Analytical, Phasor and Diffusion simulations.
%
%   Uses parameterized class test
%   for more info: 
%   https://www.mathworks.com/help/matlab/matlab_prog/create-basic-parameterized-test.html
%   
%
%========================  CORSMED AB © 2020 ==============================
%

% define the class
classdef crossComparisonTest < matlab.unittest.TestCase   
    
    % parameterization of the tests
    properties (TestParameter)
        expControl  = {data.expControl.initialize()};
        t1Value     = {1.0};
        t2Value     = {0.5};
        gradDir     = {'x', 'y', 'z'};
        gradAmp     = {0.000, 0.030}
        gradTime    = {5e-3, 10e-3};
        gradSlew    = {0, 150};
        TAU         = {25e-3, 50e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% compare all methods
        function runCrossComparison(testCase,expControl,...
                t1Value,t2Value,TAU,gradAmp,gradTime,gradSlew,gradDir)

            expControl.debug.debugMode    = 0;
            motionModel = [];
            
            %% mr System specs
            [mrSystem] = data.mrSystem.initialize();
            mrSystem.maxGStrenght   = gradAmp;
            mrSystem.SlewRate       = gradSlew;
             
            %% define RF properties
            sliceThickness  = 6e-3;
            rfDuration      = 3e-3;
            rfCycles        = 2;
            preRWScale      = 0;
            postRWScale     = 0;
            doSliceSelect   = 0;
            
            %% main RF
            mainRF      = data.acquisition.initializeRF(90);
            %% assign values
            mainRF.flipAngle          = 90; % in degrees
            mainRF.phase              = -pi/2; % in rads
            mainRF.duration           = rfDuration;
            mainRF.cycles             = rfCycles;
            mainRF.sliceThickness     = sliceThickness;
            mainRF.doSliceSelect      = doSliceSelect;
            mainRF.preRWScale         = preRWScale;
            mainRF.postRWScale        = postRWScale;
            %% refocusing RF
            refRF      = data.acquisition.initializeRF(90);
            %% assign values
            refRF.flipAngle          = 180; % in degrees
            refRF.phase              = 0; % in rads
            refRF.duration           = rfDuration;
            refRF.cycles             = rfCycles;
            refRF.sliceThickness     = sliceThickness;
            refRF.doSliceSelect      = doSliceSelect;
            refRF.preRWScale         = preRWScale;
            refRF.postRWScale        = postRWScale;
            %% PG encoding
            encPG       = data.acquisition.initializePG();
            %% assign values
            encPG.Dir = gradDir;
            encPG.TG  = gradTime;
            encPG.AG  = gradAmp;
            encPG.TAU = TAU;
            
            %% generate the sequence
            [pulseSequence] = sequence.familyPGSE.encodingPGSE( ...
                mainRF, refRF, encPG, mrSystem, expControl) ;

            % cheat: modify sequence to have 1 RO only at the end
            pulseSequence.numRxs        = 1;
            pulseSequence.rxSignal(:)   = 0;
            pulseSequence.rxSignal(end) = 1;
            pulseSequence.rxLimits      = [pulseSequence.numSteps, pulseSequence.numSteps];

            %% define conditions and initial magnetization state
            b0                  = 1.0;
            mu                  = 1;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 0.0;
            
            %% domain
            fovX = 0.010;
            fovY = 0.010;
            fovZ = 0.010;
            dx   = 1e-3;
            dy   = 1e-3;
            dz   = 1e-3;
            volFov  = fovX*fovY*fovZ;

            %% generate spin Model
            [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ,...
                dx, dy, dz, t1Value, t2Value, b0, mu, expControl.debug.debugMode);
            
            %% initial magnetizations
            % get sizes
            numCoils = spinModel.slice{1}.model.numCoils;
            numIso   = spinModel.slice{1}.model.numIsochromats;
            numRxs   = 1;
            
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
                solAnalytical, pulseSequence,spinModel.slice{1}.model,...
                motionModel, expControl);
            
            %% call the Phasor kernel
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
                solPhasor, pulseSequence,spinModel.slice{1}.model,...
                motionModel, expControl);
            
            %% call the diffusion kernel
            solDiff = data.solution.initialize(numIso,numRxs,numCoils);
            solDiff.indexes    = [];
            solDiff.Mm         = 1*M0m;
            solDiff.Mp         = 1*M0p;
            solDiff.Mx         = 1*M0x;
            solDiff.My         = 1*M0y;
            solDiff.Mz         = 1*M0z;
            solDiff.dMpDx      = 0*M0p;
            solDiff.dMpDy      = 0*M0p;
            solDiff.dMpDz      = 0*M0p;
            solDiff.Sx         = zeros(numCoils*numRxs,1);
            solDiff.Sy         = zeros(numCoils*numRxs,1);
            solDiff.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solDiff] = simulator.diffusion.kernel.runSDW(...
                solDiff, pulseSequence,spinModel.slice{1}.model,...
                motionModel, expControl);
            
            %% reference signals
            totalTime  = pulseSequence.time(end);
            refSignalX = initialMagnitude*volFov*exp(-totalTime/t2Value);
            refSignalY = 0.0;
            
            %% verify correctness with exponential
            maxErrVal  = abs(refSignalX + 1j*refSignalY)*1e-2;
            testCase.verifyLessThan(abs(solAnalytical.Sx(end)-refSignalX),...
                maxErrVal);
            testCase.verifyLessThan(abs(solAnalytical.Sy(end)-refSignalY),...
                maxErrVal);
            testCase.verifyLessThan(abs(solPhasor.Sx(end)-refSignalX),...
                maxErrVal);
            testCase.verifyLessThan(abs(solPhasor.Sy(end)-refSignalY),...
                maxErrVal);
            testCase.verifyLessThan(abs(solDiff.Sx(end)-refSignalX),...
                maxErrVal);
            testCase.verifyLessThan(abs(solDiff.Sy(end)-refSignalY),...
                maxErrVal);
            
            %% verify phasor signals with analytical
            maxErrVal  = abs(refSignalX + 1j*refSignalY)*1e-3;
            testCase.verifyLessThan(abs(solPhasor.Sx - solAnalytical.Sx),...
                maxErrVal);
            testCase.verifyLessThan(abs(solPhasor.Sy - solAnalytical.Sy),...
                maxErrVal);
            testCase.verifyLessThan(abs(solPhasor.Sz - solAnalytical.Sz),...
                maxErrVal);
            
            %% verify phasor signals with diffusion
            maxErrVal  = abs(refSignalX + 1j*refSignalY)*1e-6;
            testCase.verifyLessThan(abs(solPhasor.Sx - solDiff.Sx),...
                maxErrVal);
            testCase.verifyLessThan(abs(solPhasor.Sy - solDiff.Sy),...
                maxErrVal);
            testCase.verifyLessThan(abs(solPhasor.Sz - solDiff.Sz),...
                maxErrVal);
            
        end
        
    end
end