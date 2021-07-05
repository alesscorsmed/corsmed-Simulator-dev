%
% SIMULATOR.BLOCH.TEST.SIGNALINTEGRATIONTEST
%
%	Functional unit testing for signal integration in Bloch Simulator.
%   Verifies that the signal integration works when compared with a finely
%   discretized grid.
%
%   Uses parameterized class test
%   for more info: 
%   https://www.mathworks.com/help/matlab/matlab_prog/create-basic-parameterized-test.html
%
% Tests description
%   
%
%========================  CORSMED AB © 2020 ==============================
%

% define the class
classdef signalIntegrationTest < matlab.unittest.TestCase   
    
    % parameterization of the tests
    properties (TestParameter)
        expControl  = {data.expControl.initialize()};
        t1Value     = {0.5};
        t2Value     = {0.1, 0.25};
        gradDir     = {'x', 'y', 'z'};
        gradAmp     = {0.01, 0.03}
        gradTime    = {3e-3};
        gradSlew    = {0, 150};
        TAU         = {5e-3, 10e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        % test the signal integration
        function runSignalIntegration(testCase,expControl,...
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
            
            %% define conditions and initial magnetization state
            b0                  = 1.0;
            mu                  = 1;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 90.0;
            
            %% domain
            fovX    = 0.005;
            fovY    = 0.005;
            fovZ    = 0.005;
            volFov  = fovX*fovY*fovZ;
        
            %% reference value at end
            totalTime  = pulseSequence.time(end);
            iniSignalX = initialMagnitude*volFov*cos(initialPhase*pi/180);
            endSignalX = iniSignalX*exp(-totalTime/t2Value);
            iniSignalY = initialMagnitude*volFov*sin(initialPhase*pi/180);
            endSignalY = iniSignalY*exp(-totalTime/t2Value);
            endSignalZ = initialMagnitude*volFov*(1-exp(-totalTime/t1Value));
            
            %% generate spin Model for reference resolution
            nVox = 50;
            dx = fovX/nVox;
            dy = fovY/nVox;
            dz = fovZ/nVox;
            %% generate spin Model
            [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ,...
                dx, dy, dz, t1Value, t2Value, b0, mu, expControl.debug.debugMode);
            
            % shift the coordinates of the model for more impact of the gradient
            spinModel.slice{1}.model.r3D = spinModel.slice{1}.model.r3D + fovX;
            spinModel.slice{1}.model.x = spinModel.slice{1}.model.x + fovX;
            spinModel.slice{1}.model.y = spinModel.slice{1}.model.y + fovX;
            spinModel.slice{1}.model.z = spinModel.slice{1}.model.z + fovX;
            
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
                solAnalytical, pulseSequence,spinModel.slice{1}.model,...
                motionModel, expControl);
           
            %% generate spin Model with 1 voxel
            nVox = 1;
            dx = fovX/nVox;
            dy = fovY/nVox;
            dz = fovZ/nVox;
            [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ,...
                dx, dy, dz, t1Value, t2Value, b0, mu, expControl.debug.debugMode);
            
            % shift the coordinates of the model for more impact of the gradient
            spinModel.slice{1}.model.r3D = spinModel.slice{1}.model.r3D + fovX;
            spinModel.slice{1}.model.x = spinModel.slice{1}.model.x + fovX;
            spinModel.slice{1}.model.y = spinModel.slice{1}.model.y + fovX;
            spinModel.slice{1}.model.z = spinModel.slice{1}.model.z + fovX;
            
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
            
            %% get errors
            errorSignalX = abs(solPhasor.Sx - solAnalytical.Sx);
            errorSignalY = abs(solPhasor.Sy - solAnalytical.Sy);
            errorSignalZ = abs(solPhasor.Sz - solAnalytical.Sz);
            
            %% verify correctness with exponential
            testCase.verifyLessThan(abs(solPhasor.Sx(end)-endSignalX),...
                abs(iniSignalX + 1j*iniSignalY)*1e-2);
            testCase.verifyLessThan(abs(solPhasor.Sy(end)-endSignalY),...
                abs(iniSignalX + 1j*iniSignalY)*1e-2);
            testCase.verifyLessThan(abs(solPhasor.Sz(end)-endSignalZ),...
                abs(iniSignalX + 1j*iniSignalY)*1e-2);
            
            %% verify signals with reference analytical
            testCase.verifyLessThan(errorSignalX, abs(iniSignalX + 1j*iniSignalY)*1e-2);
            testCase.verifyLessThan(errorSignalY, abs(iniSignalX + 1j*iniSignalY)*1e-2);
            testCase.verifyLessThan(errorSignalZ, abs(iniSignalX + 1j*iniSignalY)*1e-2);

        end
    end
end