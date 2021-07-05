%
% SIMULATOR.BLOCH.TEST.BLOCHREFOCUSINGTEST
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
classdef blochRefocusingTest < matlab.unittest.TestCase   
    
    % parameterization of the tests
    properties (TestParameter)
        expControl  = {data.experiment.initializeExpControl()};
        t1Value     = {1.0, 5.0};
        t2Value     = {0.5, 2.5};
        gradDir     = {'x', 'y', 'z'};
        gradAmp     = {0.030}
        gradTime    = {5e-3, 10e-3};
        gradSlew    = {0, 150};
        TAU         = {25e-3, 50e-3};
    end
    
    % Functions with unit tests
    methods (Test)
        
        % test the the flip angle without slice selection
        function runRefocusingTest(testCase,expControl,...
                t1Value,t2Value,...
                TAU,gradAmp,gradTime,gradSlew,gradDir)
            
            expControl.debug.debugMode = 1;
            expControl.sequence.dtRF   = 1e-6;
            expControl.application = 'test';
            mrSystem = expControl.mrSystem;
            motionModel = [];
            
            %% define RF properties
            sliceThickness  = 6e-3;
            rfDuration      = 3e-3;
            rfCycles        = 2;
            preRWScale      = 0;
            postRWScale     = 0;
            doSliceSelect   = 0;
            
            %% main RF
            mainRF      = data.experiment.initializeAcqRF(90);
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
            refRF      = data.experiment.initializeAcqRF(180);
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
            encPG       = data.experiment.initializeAcqPG();
            %% assign values
            encPG.Dir = gradDir;
            encPG.TG  = gradTime;
            encPG.AG  = gradAmp;
            encPG.TAU = TAU;
            
            %% generate the sequence
            [pulseSequence] = sequence.familyPGSE.encodingPGSE( ...
                mainRF, refRF, encPG, mrSystem, expControl) ;
            
            %% define conditions and initial magnetization state
            b0                  = 1.0;
            mu                  = 1;
            initialMagnitude    = b0*mu;
            initialPhase        = 0.0;
            initialFA           = 0.0;
            
            %% domain
            fovX = 0.001;
            fovY = 0.001;
            fovZ = 0.001;
            dx   = 1e-3;
            dy   = 1e-3;
            dz   = 1e-3;
            switch lower(gradDir)
                case 'x'
                    fovX = 0.100;
                    dx = 1e-4;
                case 'y'
                    fovY = 0.100;
                    dy = 1e-4;
                case 'z'
                    fovZ = 0.100;
                    dz = 1e-4;
            end
            
            %% generate spin Model
            [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ, dx, dy, dz,...
                t1Value, t2Value, b0, mu, expControl.debug.debugMode);
            
            % get sizes
            numCoils = spinModel.slice{1}.model.numRxCoils;
            numIso   = spinModel.slice{1}.model.numIsochromats;
            numRxs   = pulseSequence.numRxs;
            
            % modifty model to have 2 tissues
            idxTissue1 = 1:round(numIso/2);
            idxTissue2 = 1+round(numIso/2):numIso;
            spinModel.slice{1}.model.numTissues         = 2;
            spinModel.slice{1}.model.numProperties      = 6;
            spinModel.slice{1}.model.tissueValues       = zeros(2,6);
            spinModel.slice{1}.model.tissueValues(1,1)  = t1Value;
            spinModel.slice{1}.model.tissueValues(1,2)  = t2Value;
            spinModel.slice{1}.model.tissueValues(2,1)  = t1Value/2;
            spinModel.slice{1}.model.tissueValues(2,2)  = t2Value/2;
            spinModel.slice{1}.model.tissueType(idxTissue1) = 1;
            spinModel.slice{1}.model.tissueType(idxTissue2) = 2;
            % tissue diffusion
            spinModel.slice{1}.model.tissueDiff = zeros(2,3);   
            
            %% initial magnetizations
            M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
            M0p = initialPhase*ones(numIso,1);
            M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
            M0x = M0m.*cos(initialPhase);
            M0y = M0m.*sin(initialPhase);
                       
            %% call the standard analytical kernel            
            solAnalytical = data.simulation.initializeSolution(numIso,numRxs,numCoils);
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
            solPhasor = data.simulation.initializeSolution(numIso,numRxs,numCoils);
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
            
            %% target magnetizations
            totalTime  = pulseSequence.time(end);
            refMm = zeros(numIso,1);
            refMm(idxTissue1) = initialMagnitude*exp(-totalTime/t2Value);
            refMm(idxTissue2) = initialMagnitude*exp(-totalTime/(t2Value/2));
            refMp = zeros(numIso,1); % by pulse construction
            
            %% verify correctness
            testCase.verifyLessThan(abs(solAnalytical.Mm-refMm),...
                initialMagnitude*1e-2);
            testCase.verifyLessThan(abs(solAnalytical.Mp-refMp),...
                1e-1);
            testCase.verifyLessThan(abs(solPhasor.Mm-refMm),...
                initialMagnitude*1e-2);
            testCase.verifyLessThan(abs(solPhasor.Mp-refMp),...
                1e-1);
            
            % verify derivatives
            % make sure the value that will be used 
            % for phase integration is close to 0
            switch lower(gradDir)
                case 'x'
                    intPhase = abs(solPhasor.dMpDx)*dx*0.5;
                    testCase.verifyLessThan(intPhase, 0.01);
                case 'y'
                    intPhase = abs(solPhasor.dMpDy)*dy*0.5;
                    testCase.verifyLessThan(intPhase, 0.01);
                case 'z'
                    intPhase = abs(solPhasor.dMpDz)*dz*0.5;
                    testCase.verifyLessThan(intPhase, 0.01);
            end
            
        end
    end
end