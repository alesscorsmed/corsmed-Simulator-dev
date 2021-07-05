%
% SPINTWIN.TEST.UNIT.FIDTEST
%
%	Functional unit testing for the Bloch Simulator.
%   Verifies correctness exponential decays
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
classdef fidTest < matlab.unittest.TestCase
    
    % parameterization of the tests
    properties (TestParameter)
        precision       = {'double'};
        odeMethod       = {'analytical', 'explicit'};
        engine          = {'bloch', 'phasor', 'diffusion'};
        partType        = {'RF', 'GR', 'RO'};
        gzAmplitude     = {0.010};
        initialPhase    = {0.0, 90.0}
    end
    
    % Functions with unit tests
    methods (Test)
        
        %% test the the FID decay
        function runFidTest(testCase,...
                precision, odeMethod, engine, partType, ...
                gzAmplitude, initialPhase )
            
            %% general info for the test
            % 10ms relaxation at 0.1us discretization
            totalime    = 10e-3;
            dt          = 1e-7;
            numSteps    = round(totalime/dt);
            t2Limits    = [0.010,  5.000];
            t1Limits    = [0.040,  5.000];
            numT2       = 128;
            numT1       = 64;
            %% initial magnetization state
            b0          = 1.0;
            mu          = 1.0;
            initialM    = b0*mu;
            initialP    = initialPhase*pi/180;
            initialFA   = 90.0; % in degrees
            
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
            
            %% empty sequence
            pulseSequence = data.simulation.initializeSequence();
            pulseSequence.type      = 'FID';
            % waveforms
            pulseSequence.time      = zeros(numSteps,1);
            pulseSequence.timeDiff  = zeros(numSteps,1);
            pulseSequence.rxSignal  = zeros(numSteps,1); % receiver readout
            pulseSequence.swcSignal = zeros(numSteps,1); % software crusher
            pulseSequence.gxSignal  = zeros(numSteps,1); % x gradient
            pulseSequence.gySignal  = zeros(numSteps,1); % y gradient
            pulseSequence.gzSignal  = zeros(numSteps,1); % z gradient
            pulseSequence.rfmSignal = zeros(numSteps,1); % RF magnitude
            pulseSequence.rfpSignal = zeros(numSteps,1); % RF phase
            pulseSequence.rffSignal = zeros(numSteps,1); % RF frequency
            pulseSequence.gdwSignal = zeros(numSteps,3); % diffusion gradients
            % assign blocks
            pulseSequence.time(:)       = dt*(1:numSteps);
            pulseSequence.timeDiff(:)   = dt;
            pulseSequence.gzSignal(:)   = gzAmplitude;
            pulseSequence.rxSignal(end) = 1;
            pulseSequence.totalTime     = dt*numSteps; % total time in seconds
            pulseSequence.numSteps      = numSteps;   % number of time steps
            pulseSequence.numRxs        = 1;   % number of readout points
            % for splitting into parts (RF vs Non-RF)
            pulseSequence.numParts      = 1; % number of parts
            pulseSequence.partType{1}   = partType; % type of part: RF / RO / GR / DF
            pulseSequence.partLimits    = [1, numSteps]; % index start/end of part
            % for indicating RO parts
            pulseSequence.rxLimits      = [numSteps, numSteps];
            
            %% model
            dx = 1e-3; dy = 1e-3; dz = 1e-3;
            fovX = numT2*dx;
            fovY = numT1*dy;
            fovZ = 10.5*dz;
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
            
            %% get dimensions and initialize magnetizations
            numCoils = simModel.numRxCoils;
            numIso   = simModel.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            gamma    = pulseSequence.gamma;
            % magnetizations
            M0m = initialM*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialP*ones(numIso,1);
            M0z = initialM*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialP);
            M0y = M0m.*sin(initialP);
            
            %% references: target result
            % cumulated area
            theta = sum(pulseSequence.timeDiff(:).*pulseSequence.gzSignal(:));
            R  = -2*pi*gamma*theta*simModel.z; % phase for each voxel
            E2 = exp(-pulseSequence.totalTime./simModel.tissueValues(simModel.tissueType(:),2));
            E1 = exp(-pulseSequence.totalTime./simModel.tissueValues(simModel.tissueType(:),1));
            M0m   = abs(M0x + 1j*M0y);
            refMx = M0m.*E2.*cos(initialP + R);
            refMy = M0m.*E2.*sin(initialP + R);
            refMm = abs(refMx + 1j*refMy);
            refMp = angle(refMx + 1j*refMy);
            refMz = initialM + ( M0z - initialM ).*E1;
            
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
            RELTOL = 1e-4;
            testCase.verifyLessThan(abs(solSim.Mx-refMx),...
                initialM*RELTOL);
            testCase.verifyLessThan(abs(solSim.My-refMy),...
                initialM*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mz-refMz),...
                initialM*RELTOL);
            testCase.verifyLessThan(abs(solSim.Mm-refMm),...
                initialM*RELTOL);
            
        end
    end
end