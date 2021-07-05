%
% SPINTWIN.TEST.UNIT.FLIPANGLETEST
%
%	Functional unit testing for the Bloch Simulator.
%   Verifies correctness of RF flipping with no slice selection
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
classdef flipAngleTest < matlab.unittest.TestCase
    
    % parameterization of the tests
    properties (TestParameter)
        precision       = {'double'};
        odeMethod       = {'analytical', 'explicit'};
        engine          = {'bloch', 'phasor', 'diffusion'};
        rfDuration      = {1e-3, 3e-3};
        rfFlipAngle     = {45, 90, 180};
        rfPhase         = {0.0, 90, -90, 180};
        rfCycles        = {2, 4};  
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test the the flip angle without slice selection
        % Compare with theoretical result (with loose accuracy)
        function theoreticalFlipAngle(testCase,...
                precision, odeMethod, engine, ...
                rfDuration, rfFlipAngle, rfPhase, rfCycles )
        
            %% prepare the data required
            % mrSystem
            mrSystem.b0             = 1.5000;
            mrSystem.maxGStrenght   = 0.0400;
            mrSystem.SlewRate       = 0.0;
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
            % RF data
            rfData.type             = 'sinc';
            rfData.flipAngle      	= rfFlipAngle; % in degrees
            rfData.phase         	= rfPhase*pi/180; % in rads
            rfData.duration         = rfDuration;
            rfData.cycles           = rfCycles;
            rfData.sliceThickness   = 6e-3;
            rfData.doSliceSelect 	= 0;
            rfData.preRWScale    	= 0;
            rfData.postRWScale    	= 0;
            rfData.keepTimingSS     = 0;
            % sequence
            [ pulseSequence ] = spinTwin.test.seq.singleRF( ...
                rfData, mrSystem, dt);
            
            %% model
            b0 = 1.0;
            mu = 1.0;
            dx = 1e-3;
            dy = 1e-3;
            dz = 1e-4;
            nZ = 2*(rfData.sliceThickness/dz)+1; % make sure is odd to have 0 position
            fovX = dx;
            fovY = dy;
            fovZ = nZ*dz;
            t1Value = 1e3;
            t2Value = 5e3;
            % generate spin Model
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            % get and modify the simulation model
            simModel = spinModel.slice{1}.model;
            
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

            %% target (ideal) magnetizations
            refMm = initialMagnitude*sin(rfFlipAngle*pi/180)*ones(numIso,1);
            refMp = initialPhase*ones(numIso,1) + rfPhase*pi/180 + pi/2;
            refMz = initialMagnitude*cos(rfFlipAngle*pi/180)*ones(numIso,1);
            refMx = refMm.*cos(refMp);
            refMy = refMm.*sin(refMp);
        
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
            
            %% verify correctness with exponential
            RELTOL = 1e-3;
            testCase.verifyLessThan(abs(solSim.Mx-refMx),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.My-refMy),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mz-refMz),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mm-refMm),...
                initialMagnitude*RELTOL);
            
        end
        
        %% test the the flip angle without slice selection
        % Compare with reference result (fine time resolution)
        function referenceFlipAngle(testCase,...
                precision, odeMethod, engine, ...
                rfDuration, rfFlipAngle, rfPhase, rfCycles )
        
            %% prepare the data required
            % mrSystem
            mrSystem.b0             = 1.5000;
            mrSystem.maxGStrenght   = 0.0400;
            mrSystem.SlewRate       = 0.0;
            % motion
            motionModel             = [];
            % dbg
            dbgControl.mode         = 0;
            dbgControl.file         = [];
            dbgControl.showWaitBar  = 0;
            %% initialize the sim Control with the default versioning for testing
            [simControl] = spinTwin.setup.initializeSimControl();
        
            %% get the sequence
            % RF data
            rfData.type             = 'sinc';
            rfData.flipAngle      	= rfFlipAngle; % in degrees
            rfData.phase         	= rfPhase*pi/180; % in rads
            rfData.duration         = rfDuration;
            rfData.cycles           = rfCycles;
            rfData.sliceThickness   = 6e-3;
            rfData.doSliceSelect 	= 0;
            rfData.preRWScale    	= 0;
            rfData.postRWScale    	= 0;
            rfData.keepTimingSS     = 0;
            
            %% model
            b0 = 1.0;
            mu = 1.0;
            % T1/T2 different values
            t2Limits = [0.010,  5.000];
            t1Limits = [0.040,  5.000];
            numT2    = 32;
            numT1    = 32;
            dx = 1e-3; dy = 1e-3; dz = 1e-3;
            nZ = 1;
            fovX = numT2*dx;
            fovY = numT1*dy;
            fovZ = nZ*dz;
            t1Value = mean(t1Limits);
            t2Value = mean(t2Limits);
            % generate spin Model
            [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
                fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
            % get and modify the simulation model
            simModel = spinModel.slice{1}.model;
            [nX, nY, nZ, ~] = size(simModel.r3D);
            numTissues = nX*nY*nZ;
            % modify the T1/T2 values to be different for each voxel
            t2Values = logspace(log10(min(t2Limits)), log10(max(t2Limits)), nX);
            t1Values = logspace(log10(min(t1Limits)), log10(max(t1Limits)), nY);
            % nX x nY with different T1/T2
            t2Mat = repmat(reshape(t2Values,nX,1),1,nY);
            t1Mat = t2Mat + repmat(reshape(t1Values,1,nY),nX,1); % make sure T1 > T2
            % extend to nZ
            t2Mat = repmat(t2Mat,1,nZ);
            t1Mat = repmat(t1Mat,1,nZ);
            tissueType = reshape(1:numTissues,[],1);
            % assign
            simModel.numTissues         = numTissues;
            simModel.numProperties      = 6;
            simModel.tissueValues       = zeros(numTissues,6);
            simModel.tissueValues(:,1)  = t1Mat(:);
            simModel.tissueValues(:,2)  = t2Mat(:);
            simModel.tissueType(:)      = reshape(tissueType,[],1);
            % assign zero diffusion
            simModel.tissueDiff         = zeros(numTissues,3);

            %% reference sequence: dt = 1e-8
            % get the sequence
            dt = 1e-8;
            [ pulseSequenceRef ] = spinTwin.test.seq.singleRF( ...
                rfData, mrSystem, dt);
            
             %% testing sequence: dt = 1e-7
             dt = 1e-7;
            [ pulseSequence ] = spinTwin.test.seq.singleRF( ...
                rfData, mrSystem, dt);
                        
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
            
            %% reference run: double precision at 1e-7 (ref sequence)
            % call the kernel
            simControl.precision        = 'double'; % single / double
            simControl.simulationEngine = 'Bloch';
            simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
            % initialize solution
            solRef = data.simulation.initializeSolution(numIso,numRxs,numCoils);
            solRef.indexes    = [];
            solRef.Mx         = 1*M0x;
            solRef.My         = 1*M0y;
            solRef.Mz         = 1*M0z;
            solRef.dMpDx      = 0*M0x;
            solRef.dMpDy      = 0*M0x;
            solRef.dMpDz      = 0*M0x;
            solRef.Sx         = zeros(numCoils*numRxs,1);
            solRef.Sy         = zeros(numCoils*numRxs,1);
            solRef.Sz         = zeros(numCoils*numRxs,1);
            % kernel call
            [solRef, ~] = spinTwin.fwdBloch.runBloch(...
                solRef, pulseSequenceRef, simModel,...
                motionModel, simControl, dbgControl );
            % magnitude and phase
            solRef.Mm = abs(solRef.Mx + 1j*solRef.My);
            solRef.Mp = angle(solRef.Mx + 1j*solRef.My);
            
            %% call the testing kernel
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
            
            %% verify correctness with exponential
            RELTOL = 1e-4;
            testCase.verifyLessThan(abs(solSim.Mx-solRef.Mx),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.My-solRef.My),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mz-solRef.Mz),...
                initialMagnitude*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mm-solRef.Mm),...
                initialMagnitude*RELTOL);
            
        end      
        
    end
end