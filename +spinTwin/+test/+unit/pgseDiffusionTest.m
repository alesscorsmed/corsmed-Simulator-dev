%
% SPINTWIN.TEST.UNIT.PGSEDIFFUSIONTEST
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
        precision       = {'double'};
        odeMethod       = {'analytical', 'explicit'};
        engine          = {'diffusion'};
        gradDir         = {'x', 'y', 'z'};
        gradAmp         = {0.001, 0.030};
        gradSlew        = {0.0, 150.0};
        gradDuration    = {5e-3, 10e-3, 20e-3};
        echoTime        = {25e-3, 50e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test with no RF (reversed gradients)
        function runPGReversed(testCase, ...
                precision, engine, ...
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

            %% domain and model
            b0 	 = 1.0;
            mu 	 = 1.0;
            fovX = 0.010;
            fovY = 0.010;
            fovZ = 0.010;
            dx   = 5e-4;
            dy   = 5e-4;
            dz   = 5e-4;
            t1Value = 10.500;
            t2Value = 10.300;
            % generate spin Model
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            
            simModel = spinModel.slice{1}.model;
            %% get dimensions and initialize magnetizations
            numCoils = simModel.numRxCoils;
            numIso   = simModel.numIsochromats;
            
            %% modify the diffusion values to be random for each voxel
            % within limits 1e-9 to 7e-9
            xDiffusion = (1 + 6*rand(numIso,1))*1e-9;
            yDiffusion = (1 + 6*rand(numIso,1))*1e-9;
            zDiffusion = (1 + 6*rand(numIso,1))*1e-9;
            % tissue type per voxel
            tissueType = reshape(1:numIso,[],1);
            % assign
            simModel.numTissues         = numIso;
            simModel.numProperties      = 6;
            simModel.tissueValues       = zeros(numIso,6);
            simModel.tissueValues(:,1)  = t1Value;
            simModel.tissueValues(:,2)  = t2Value;
            simModel.tissueType(:)      = reshape(tissueType,[],1);
            % assign diffusion
            simModel.tissueDiff      = zeros(numIso,3);
            simModel.tissueDiff(:,1) = xDiffusion(:);
            simModel.tissueDiff(:,2) = yDiffusion(:);
            simModel.tissueDiff(:,3) = zDiffusion(:);
            
            %% initial magnetizations
            initialMagnitude = b0*mu;
            initialPhase 	 = 0.0;
            initialFA      	 = 90.0;
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
            
            %% get the sequence
            [pulseSequence, grData] = spinTwin.test.seq.pgReversed( ...
                grData, mrSystem, dt);
            
            numRxs   = pulseSequence.numRxs;
            
            %% compute the reference ADC
            beta = grData.beta;
            switch lower(gradDir)
                case 'x'
                    Dcoeff = simModel.tissueDiff(:,1);
                case 'y'
                    Dcoeff = simModel.tissueDiff(:,2);
                case 'z'
                    Dcoeff = simModel.tissueDiff(:,3);
            end
            ADCref = exp( -beta*abs(Dcoeff) );
            
            %% call the standard analytical kernel
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
                solBloch, pulseSequence,simModel,...
                motionModel, simControl, dbgControl );
            % magnitude and phase
            solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
            solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);
            
            %% call the Diffusion kernel (no RF, only w/ Analytical)
            simControl.simulationEngine = engine;
            simControl.odeMethod        = 'analytical';
            simControl.precision        = precision;
            % initialize solution
            solDiffusion = data.simulation.initializeSolution(numIso,numRxs,numCoils);
            solDiffusion.indexes    = [];
            solDiffusion.Mx         = 1*M0x;
            solDiffusion.My         = 1*M0y;
            solDiffusion.Mz         = 1*M0z;
            solDiffusion.dMpDx      = 0*M0x;
            solDiffusion.dMpDy      = 0*M0x;
            solDiffusion.dMpDz      = 0*M0x;
            solDiffusion.Sx         = zeros(numCoils*numRxs,1);
            solDiffusion.Sy         = zeros(numCoils*numRxs,1);
            solDiffusion.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solDiffusion] = spinTwin.fwdBloch.runDiffusion(...
                solDiffusion, pulseSequence,simModel,...
                motionModel, simControl, dbgControl );
            % magnitude and phase
            solDiffusion.Mm = abs(solDiffusion.Mx + 1j*solDiffusion.My);
            solDiffusion.Mp = angle(solDiffusion.Mx + 1j*solDiffusion.My);
            
            %% compute numerical ADC
            ADCsim = solDiffusion.Mm./solBloch.Mm;
            
            %% error
            ADCerror = abs(ADCsim-ADCref);
            
            %% verify signals with reference analytical
            testCase.verifyLessThan(ADCerror, 1e-2);
            
        end
        
        
        %% test with RF
        function runPGSE(testCase, ...
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
            dt                      = 1e-6;
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
            
            %% domain and model
            b0 	 = 1.0;
            mu 	 = 1.0;
            fovX = 0.010;
            fovY = 0.010;
            fovZ = 0.010;
            dx   = 2e-4;
            dy   = 2e-4;
            dz   = 2e-4;
            t1Value = 10.500;
            t2Value = 10.300;
            % generate spin Model
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            
            simModel = spinModel.slice{1}.model;
            %% get dimensions and initialize magnetizations
            numCoils = simModel.numRxCoils;
            numIso   = simModel.numIsochromats;
            
            %% modify the diffusion values to be random for each voxel
            % within limits 1e-9 to 7e-9
            xDiffusion = (1 + 6*rand(numIso,1))*1e-9;
            yDiffusion = (1 + 6*rand(numIso,1))*1e-9;
            zDiffusion = (1 + 6*rand(numIso,1))*1e-9;
            % tissue type per voxel
            tissueType = reshape(1:numIso,[],1);
            % assign
            simModel.numTissues         = numIso;
            simModel.numProperties      = 6;
            simModel.tissueValues       = zeros(numIso,6);
            simModel.tissueValues(:,1)  = t1Value;
            simModel.tissueValues(:,2)  = t2Value;
            simModel.tissueType(:)      = reshape(tissueType,[],1);
            % assign diffusion
            simModel.tissueDiff      = zeros(numIso,3);
            simModel.tissueDiff(:,1) = xDiffusion(:);
            simModel.tissueDiff(:,2) = yDiffusion(:);
            simModel.tissueDiff(:,3) = zDiffusion(:);
            
            %% initial magnetizations
            initialMagnitude = b0*mu;
            initialPhase 	 = 0.0;
            initialFA      	 = 0.0;
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
            
            %% get the sequence
            [ pulseSequence, grData ] = spinTwin.test.seq.pgse( ...
                grData, rfData, mrSystem, dt);
            % readouts
            numRxs   = pulseSequence.numRxs;
            
            %% compute the reference ADC
            beta = grData.beta;
            switch lower(gradDir)
                case 'x'
                    Dcoeff = simModel.tissueDiff(:,1);
                case 'y'
                    Dcoeff = simModel.tissueDiff(:,2);
                case 'z'
                    Dcoeff = simModel.tissueDiff(:,3);
            end
            ADCref = exp( -beta*abs(Dcoeff) );
            
            %% call the standard analytical kernel
            simControl.simulationEngine = 'Bloch';
            simControl.odeMethod        = odeMethod; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
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
                solBloch, pulseSequence,simModel,...
                motionModel, simControl, dbgControl );
            % magnitude and phase
            solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
            solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);
            
            %% call the Diffusion kernel 
            simControl.simulationEngine = engine;
            simControl.odeMethod        = odeMethod;
            simControl.precision        = precision;
            % initialize solution
            solDiffusion = data.simulation.initializeSolution(numIso,numRxs,numCoils);
            solDiffusion.indexes    = [];
            solDiffusion.Mx         = 1*M0x;
            solDiffusion.My         = 1*M0y;
            solDiffusion.Mz         = 1*M0z;
            solDiffusion.dMpDx      = 0*M0x;
            solDiffusion.dMpDy      = 0*M0x;
            solDiffusion.dMpDz      = 0*M0x;
            solDiffusion.Sx         = zeros(numCoils*numRxs,1);
            solDiffusion.Sy         = zeros(numCoils*numRxs,1);
            solDiffusion.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solDiffusion] = spinTwin.fwdBloch.runDiffusion(...
                solDiffusion, pulseSequence,simModel,...
                motionModel, simControl, dbgControl );
            % magnitude and phase
            solDiffusion.Mm = abs(solDiffusion.Mx + 1j*solDiffusion.My);
            solDiffusion.Mp = angle(solDiffusion.Mx + 1j*solDiffusion.My);
            
            %% compute numerical ADC
            ADCsim = solDiffusion.Mm./solBloch.Mm;
            
            %% error
            ADCerror = abs(ADCsim-ADCref);
            
            %% verify signals with reference analytical
            testCase.verifyLessThan(ADCerror, 5e-2);
            
        end   
            
    end
end
