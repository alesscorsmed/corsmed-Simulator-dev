%
% SPINTWIN.TEST.UNIT.SIGNALINTEGRATIONTEST
%
%	Functional unit testing for signal integration in Bloch Simulator.
%   Verifies that the signal integration works when compared with a finely
%   discretized grid.
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
classdef signalIntegrationTest < matlab.unittest.TestCase
    
    % parameterization of the tests
    properties (TestParameter)
        precision       = {'double'};
        odeMethod       = {'analytical'};
        engine          = {'phasor', 'diffusion'};
        t1Value     	= {1.000};
        t2Value         = {0.100};
        gradDir         = {'z'};
        gradAmp         = {0.010, 0.020};
        gradSlew        = {0.0, 150.0};
        gradDuration    = {3e-3, 5e-3};
        echoTime        = {15e-3, 25e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test that the signal integration  
        function runSignalIntegration(testCase,...
                precision, odeMethod, engine, ...
                t1Value, t2Value, ...
                gradDir, gradAmp, gradSlew, ....
                gradDuration, echoTime )
            
            %% prepare the data required
            % mrSystem
            mrSystem.b0             = 1.5000;
            mrSystem.maxGStrenght   = 0.040;
            mrSystem.SlewRate       = gradSlew;
            % motion
            motionModel             = [];
            % dbg
            dbgControl.mode         = 0;
            dbgControl.file         = [];
            dbgControl.showWaitBar  = 0;
            %% initialize the sim Control with the default versioning for testing
            [simControl] = spinTwin.setup.initializeSimControl();
            
            %% get the sequence
            dt                      = 1e-6;
            % gradient and echo time
            grData.echoTime         = echoTime;
            grData.duration      	= gradDuration;
            grData.amplitude     	= gradAmp;
            grData.direction        = gradDir; % encoding direction 'x' / 'y' / 'z'           
            %% get the sequence
            [pulseSequence, grData] = spinTwin.test.seq.pgReversed( ...
                grData, mrSystem, dt);
            
            %% domain
            fovX    = 0.005;
            fovY    = 0.005;
            fovZ    = 0.005;
            volFov  = fovX*fovY*fovZ;
            dx      = 1e-3;
            dy      = 1e-3;
            dz      = 1e-3;
            
            %% initial magnetizations
            b0                  = 1.0;
            mu                  = 1.0;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 90.0;
            
            %% reference value for signal at end fo simulation
            totalTime  = pulseSequence.time(end);
            iniSignalX = initialMagnitude*volFov*cos(initialPhase*pi/180);
            endSignalX = iniSignalX*exp(-totalTime/t2Value);
            iniSignalY = initialMagnitude*volFov*sin(initialPhase*pi/180);
            endSignalY = iniSignalY*exp(-totalTime/t2Value);
            endSignalZ = initialMagnitude*volFov*(1-exp(-totalTime/t1Value));

            
            %% generate spin Model for fine resolution (50 voxels)
            % discretize on the gradient direction
            nVox = 5000;
            switch lower(grData.direction)
                case 'x'
                    dx = fovX/nVox;
                case 'y'
                    dy = fovY/nVox;
                case 'z'
                    dz = fovZ/nVox;
            end
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            % get and modify the model
            simModel        = spinModel.slice{1}.model;
            % shift the coordinates of the model for more impact of the gradient
            simModel.r3D    = simModel.r3D + fovX;
            simModel.x      = simModel.x + fovX;
            simModel.y      = simModel.y + fovX;
            simModel.z      = simModel.z + fovX;
            
            %% get dimensions and initialize magnetizations
            numCoils = simModel.numRxCoils;
            numIso   = simModel.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            
            %% initial magnetizations
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
            
            %% call the standard Bloch kernel, no signal integration
            solBloch = data.simulation.initializeSolution(numIso,numRxs,numCoils);
            solBloch.indexes    = [];
            solBloch.Mx         = 1*M0x;
            solBloch.My         = 1*M0y;
            solBloch.Mz         = 1*M0z;
            solBloch.dMpDx      = 0*M0p;
            solBloch.dMpDy      = 0*M0p;
            solBloch.dMpDz      = 0*M0p;
            solBloch.Sx         = zeros(numCoils*numRxs,1);
            solBloch.Sy         = zeros(numCoils*numRxs,1);
            solBloch.Sz         = zeros(numCoils*numRxs,1);
            % standard w/ no signal integration
            simControl.simulationEngine = 'Bloch';
            % kernel call
            [solBloch] = spinTwin.fwdBloch.runBloch(...
                solBloch, pulseSequence,simModel,...
                motionModel, simControl, dbgControl );

            %% generate spin Model with 1 voxel for signal integration
            nVox = 1;
            dx = fovX/nVox;
            dy = fovY/nVox;
            dz = fovZ/nVox;
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            % get and modify the model
            simModel        = spinModel.slice{1}.model;
            % shift the coordinates of the model for more impact of the gradient
            simModel.r3D    = simModel.r3D + fovX;
            simModel.x      = simModel.x + fovX;
            simModel.y      = simModel.y + fovX;
            simModel.z      = simModel.z + fovX;
            
            %% get dimensions and initialize magnetizations
            numCoils = simModel.numRxCoils;
            numIso   = simModel.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            
            %% initial magnetizations
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
            
            %% call the kernel
            % update the controls
            simControl.simulationEngine = engine;
            simControl.odeMethod        = odeMethod;
            simControl.precision        = precision;
            % initialize solution
            solSim = data.simulation.initializeSolution(numIso,numRxs,numCoils);
            solSim.indexes    = [];
            solSim.Mx         = 1*M0x;
            solSim.My         = 1*M0y;
            solSim.Mz         = 1*M0z;
            solSim.dMpDx      = 0*M0x;
            solSim.dMpDy      = 0*M0x;
            solSim.dMpDz      = 0*M0x;
            solSim.Sx         = zeros(numCoils*numRxs,1);
            solSim.Sy         = zeros(numCoils*numRxs,1);
            solSim.Sz         = zeros(numCoils*numRxs,1);
            
            % kernel call
            switch lower(simControl.simulationEngine)
                case lower('phasor')
                    [solSim, ~] = spinTwin.fwdBloch.runPhasor( ...
                        solSim, pulseSequence, simModel,...
                        motionModel, simControl, dbgControl );
                case lower('diffusion')
                    [solSim, ~] = spinTwin.fwdBloch.runDiffusion(...
                        solSim, pulseSequence, simModel,...
                        motionModel, simControl, dbgControl );
                otherwise % basic Bloch
                    [solSim, ~] = spinTwin.fwdBloch.runBloch(...
                        solSim, pulseSequence, simModel,...
                        motionModel, simControl, dbgControl );
            end
            
            %% get errors
            errorSignalX = abs(solSim.Sx - solBloch.Sx);
            errorSignalY = abs(solSim.Sy - solBloch.Sy);
            
            %% verify correctness with exponential
            RELTOL = 1e-4;
            ABSTOL = initialMagnitude*volFov*RELTOL;
            testCase.verifyLessThan(abs(solSim.Sx(end)-endSignalX),...
                ABSTOL);
            testCase.verifyLessThan(abs(solSim.Sy(end)-endSignalY),...
                ABSTOL);
            testCase.verifyLessThan(abs(solSim.Sz(end)-endSignalZ),...
                ABSTOL);
            
            %% verify signals with reference analytical
            RELTOL = 1e-3;
            ABSTOL = initialMagnitude*volFov*RELTOL;
            testCase.verifyLessThan(errorSignalX, ABSTOL);
            testCase.verifyLessThan(errorSignalY, ABSTOL);

        end
    end
end