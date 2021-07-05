%
% SPINTWIN.TEST.UNIT.PHASEREFOCUSTEST
%
%	Functional unit testing for the Bloch Simulator.
%   Verifies correct 180 RF refocusing
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
classdef phaseRefocusTest < matlab.unittest.TestCase
    
    % parameterization of the tests
    properties (TestParameter)
        precision       = {'double'};
        odeMethod       = {'analytical', 'explicit'};
        engine          = {'phasor', 'diffusion'};
        gradDir         = {'x', 'y', 'z'};
        gradAmp         = {0.030};
        gradSlew        = {0.0, 150.0};
        gradDuration    = {5e-3, 10e-3};
        echoTime        = {25e-3, 50e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test that the phase derivative after a refocusing 180 is zero 
        function runPhaseRefocus(testCase,...
                precision, odeMethod, engine, ...
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
            dt                      = 1e-7;
            % gradient and echo time
            grData.echoTime         = echoTime;
            grData.duration      	= gradDuration;
            grData.amplitude     	= gradAmp;
            grData.direction        = gradDir; % encoding direction 'x' / 'y' / 'z'
            % RF data
            rfData.type             = 'sinc';
            rfData.flipAngle      	= 90; % in degrees
            rfData.phase         	= 0.0; % in rads
            rfData.duration         = 1e-3;
            rfData.cycles           = 2;
            rfData.sliceThickness   = 3e-3;
            rfData.doSliceSelect 	= 0;
            rfData.preRWScale    	= 0;
            rfData.postRWScale    	= 0;
            rfData.keepTimingSS     = 0;
            
            %% get the sequence
            [ pulseSequence, grData ] = spinTwin.test.seq.pgse( ...
                grData, rfData, mrSystem, dt);
            
            %% domain and model
            b0 	 = 1.0;
            mu 	 = 1.0;
            fovX = 1e-3;
            fovY = 1e-3;
            fovZ = 1e-3;
            dx   = 1e-3;
            dy   = 1e-3;
            dz   = 1e-3;
            t1Value = 10.500;
            t2Value = 10.300;
            % discretize on the gradient direction
            switch lower(grData.direction)
                case 'x'
                    fovX = 0.050; dx = 1e-4;
                case 'y'
                    fovY = 0.050; dy = 1e-4;
                case 'z'
                    fovZ = 0.050; dz = 1e-4;
            end
            % generate spin Model
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            % get the simulation model
            simModel = spinModel.slice{1}.model;
            % voxel positions
            switch lower(grData.direction)
                case 'x'
                    position = simModel.x;
                case 'y'
                    position = simModel.y;
                case 'z'
                    position = simModel.z;
            end
            
            %% get dimensions and initialize magnetizations
            numCoils = simModel.numRxCoils;
            numIso   = simModel.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            
            %% initial magnetizations
            initialMagnitude = b0*mu;
            initialPhase 	 = 0.0;
            initialFA      	 = 0.0;
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
            
            %% call the standard analytical kernel for reference
            simControl.simulationEngine = 'Bloch';
            simControl.odeMethod        = 'analytical'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
            simControl.precision        = 'double';
            % initialize solution
            solBloch = data.simulation.initializeSolution(numIso,numRxs,numCoils);
            solBloch.indexes    = [];
            solBloch.Mx         = 1*M0x;
            solBloch.My         = 1*M0y;
            solBloch.Mz         = 1*M0z;
            solBloch.dMpDx      = 0*M0x;
            solBloch.dMpDy      = 0*M0x;
            solBloch.dMpDz      = 0*M0x;
            solBloch.Sx         = zeros(numCoils*numRxs,1);
            solBloch.Sy         = zeros(numCoils*numRxs,1);
            solBloch.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solBloch] = spinTwin.fwdBloch.runBloch(...
                solBloch, pulseSequence,spinModel.slice{1}.model,...
                motionModel, simControl, dbgControl );
            % magnitude and phase
            solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
            solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);
            
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
            
            % magnitude and phase
            solSim.Mm = abs(solSim.Mx + 1j*solSim.My);
            solSim.Mp = angle(solSim.Mx + 1j*solSim.My);
            
            %% verify correctness with reference
            RELTOL = 1e-4;
            testCase.verifyLessThan(abs(solSim.Mx-solBloch.Mx),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.My-solBloch.My),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mz-solBloch.Mz),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mm-solBloch.Mm),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mp-solBloch.Mp),...
                RELTOL);
            
            % verify derivatives
            % make sure the value that will be used for phase integration is close to 0
            ABSTOL = 1e-4;
            switch lower(gradDir)
                case 'x'
                    intPhase = abs(solSim.dMpDx)*dx*0.5;
                    testCase.verifyLessThan(intPhase, ABSTOL);
                case 'y'
                    intPhase = abs(solSim.dMpDy)*dy*0.5;
                    testCase.verifyLessThan(intPhase, ABSTOL);
                case 'z'
                    intPhase = abs(solSim.dMpDz)*dz*0.5;
                    testCase.verifyLessThan(intPhase, ABSTOL);
            end
            
        end
        
    end
end