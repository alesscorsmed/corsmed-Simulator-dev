%
% SIMULATOR.DIFUSSION.TEST.PGSEDIFFUSIONTEST
%
%	Functional unit testing for self diffusion in Bloch Simulator.
%   Verifies that the ADC matches the theoretical.
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
classdef pgseDiffusionTest < matlab.unittest.TestCase   
    
    % parameterization of the tests
    properties (TestParameter)
        expControl  = {data.expControl.initialize()};
        t1Value     = {0.7};
        t2Value     = {0.5};
        gradDir     = {'x', 'y', 'z'};
        gradAmp     = {0.001, 0.010, 0.030}
        gradTime    = {5e-3, 10e-3, 20e-3};
        gradSlew    = {0, 150};
        TAU         = {25e-3, 50e-3, 100e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test with no RF (reversed gradients)
        function runPGReversed(testCase,expControl,...
                t1Value,t2Value,TAU,gradAmp,gradTime,gradSlew,gradDir)

            expControl.debug.debugMode    = 0;
            motionModel = [];
            
            %% mr System specs
            [mrSystem] = data.mrSystem.initialize();
            mrSystem.maxGStrenght   = gradAmp;
            mrSystem.SlewRate       = gradSlew;

            %% PG encoding
            encPG       = data.acquisition.initializePG();
            %% assign values
            encPG.Dir = gradDir;
            encPG.TG  = gradTime;
            encPG.AG  = gradAmp;
            encPG.TAU = TAU;
            
            %% populate the sequence
            [pulseSequence] = sequence.familyPGSE.reversedPG(...
                encPG, mrSystem, expControl);
            
            % cheat: modify sequence to have 1 RO only at the end
            pulseSequence.numRxs        = 1;
            pulseSequence.rxSignal(:)   = 0;
            pulseSequence.rxSignal(end) = 1;
            
            %% define conditions and initial magnetization state
            b0                  = 1.0;
            mu                  = 1;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 90.0;
            
            %% domain
            fovX = 0.010;
            fovY = 0.010;
            fovZ = 0.010;
            dx   = 5e-4;
            dy   = 5e-4;
            dz   = 5e-4;
            
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
            
            %% modify the diffusion values to be random for each voxel
            % random 
            spinModel.slice{1}.model.numTissues         = numIso;
            spinModel.slice{1}.model.numProperties      = 6;
            spinModel.slice{1}.model.tissueValues       = zeros(numIso,6);
            spinModel.slice{1}.model.tissueValues(:,1)  = t1Value;
            spinModel.slice{1}.model.tissueValues(:,2)  = t2Value;
            spinModel.slice{1}.model.tissueType         = reshape(1:numIso,numIso,1);
            % tissue diffusion
            spinModel.slice{1}.model.tissueDiff = zeros(numIso,3);            
            spinModel.slice{1}.model.tissueDiff(:,1) = (1 + 6*rand(numIso,1))*1e-9;
            spinModel.slice{1}.model.tissueDiff(:,2) = (1 + 6*rand(numIso,1))*1e-9;
            spinModel.slice{1}.model.tissueDiff(:,3) = (1 + 6*rand(numIso,1))*1e-9;
            % diffusion arrays
            spinModel.slice{1}.model.xDiffusion = ...
                spinModel.slice{1}.model.tissueDiff(spinModel.slice{1}.model.tissueType,1);
            spinModel.slice{1}.model.yDiffusion = ...
                spinModel.slice{1}.model.tissueDiff(spinModel.slice{1}.model.tissueType,2);
            spinModel.slice{1}.model.zDiffusion = ...
                spinModel.slice{1}.model.tissueDiff(spinModel.slice{1}.model.tissueType,3);
            
            %% compute the reference ADC
            beta = pulseSequence.encPG.Beta;
            switch lower(gradDir)
                case 'x'
                    Dcoeff = spinModel.slice{1}.model.xDiffusion;
                case 'y'
                    Dcoeff = spinModel.slice{1}.model.yDiffusion;
                case 'z'
                    Dcoeff = spinModel.slice{1}.model.zDiffusion;
            end
            ADCref = exp( -beta*abs(Dcoeff) );
            
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
            
            %% compute numerical ADC
            ADCsim = solDiff.Mm./solAnalytical.Mm;
            
            %% error
            ADCerror = abs(ADCsim-ADCref);
            
            %% verify signals with reference analytical
            testCase.verifyLessThan(ADCerror, 5e-2);

        end
        
        
        %% test with RF
        function runPGSE(testCase,expControl,...
                t1Value,t2Value,TAU,gradAmp,gradTime,gradSlew,gradDir)

            expControl.debug.debugMode = 0;
            expControl.sequence.dt = 1e-6;
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
            dx   = 5e-4;
            dy   = 5e-4;
            dz   = 5e-4;
            
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
            
            %% modify the diffusion values to be random for each voxel
            % random within 1e-9 to 7e-9
            spinModel.slice{1}.model.numTissues         = numIso;
            spinModel.slice{1}.model.numProperties      = 6;
            spinModel.slice{1}.model.tissueValues       = zeros(numIso,6);
            spinModel.slice{1}.model.tissueValues(:,1)  = t1Value;
            spinModel.slice{1}.model.tissueValues(:,2)  = t2Value;
            spinModel.slice{1}.model.tissueType         = reshape(1:numIso,numIso,1);
            % tissue diffusion
            spinModel.slice{1}.model.tissueDiff = zeros(numIso,3);            
            spinModel.slice{1}.model.tissueDiff(:,1) = (1 + 6*rand(numIso,1))*1e-9;
            spinModel.slice{1}.model.tissueDiff(:,2) = (1 + 6*rand(numIso,1))*1e-9;
            spinModel.slice{1}.model.tissueDiff(:,3) = (1 + 6*rand(numIso,1))*1e-9;
            % diffusion arrays
            spinModel.slice{1}.model.xDiffusion = ...
                spinModel.slice{1}.model.tissueDiff(spinModel.slice{1}.model.tissueType,1);
            spinModel.slice{1}.model.yDiffusion = ...
                spinModel.slice{1}.model.tissueDiff(spinModel.slice{1}.model.tissueType,2);
            spinModel.slice{1}.model.zDiffusion = ...
                spinModel.slice{1}.model.tissueDiff(spinModel.slice{1}.model.tissueType,3);
            
            %% compute the reference ADC
            beta = pulseSequence.encPG.Beta;
            switch lower(gradDir)
                case 'x'
                    Dcoeff = spinModel.slice{1}.model.xDiffusion;
                case 'y'
                    Dcoeff = spinModel.slice{1}.model.yDiffusion;
                case 'z'
                    Dcoeff = spinModel.slice{1}.model.zDiffusion;
            end
            ADCref = exp( -beta*abs(Dcoeff) );
            
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
            
            %% compute numerical ADC
            ADCsim = solDiff.Mm./solAnalytical.Mm;
            
            %% error
            ADCerror = abs(ADCsim-ADCref);
            
            %% verify signals with reference analytical
            testCase.verifyLessThan(ADCerror, 5e-2);

        end
               
    end
end