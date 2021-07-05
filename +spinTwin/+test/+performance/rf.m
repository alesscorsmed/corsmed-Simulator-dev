%
% SPINTWIN.TEST.PERFORMANCE.RF
%
%	Runs an RF block and test the performance an accuracy
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

% cases to run
sliceSelectList     = {'ON', 'OFF'};
precisionList       = {'double', 'single'};
% the next two go together
odeMethodList        = {'analytical', 'explicit',  'analytical', 'explicit', 'analytical', 'explicit', 'adaptiveExp', 'implicit', 'adaptiveImp', };
engineList           = {'diffusion',  'diffusion', 'phasor',     'phasor',   'bloch',      'bloch',    'bloch',       'bloch',    'bloch'};


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

%% initialize the sim Control with the default versioning for testing
[simControl] = spinTwin.setup.initializeSimControl();

%% prepare data
dt          = 1e-6;
t2Limits    = [0.010,  5.000];
t1Limits    = [0.040,  5.000];
numT2       = 256;
numT1       = 64;

%% model
b0 = 1.0;
mu = 1.0;
dx = 1e-3;
dy = 1e-3;
dz = 1e-4;
sliceThickness = 6e-3;
nZ = 2*(sliceThickness/dz)+0.5; % make sure is odd to have 0 position
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



%% loop on slice selection case
for ss = 1:length(sliceSelectList)
    
    
    %% RF
    rfData.type             = 'sinc';
    rfData.flipAngle      	= 90; % in degrees
    rfData.phase         	= 0.0; % in rads
    rfData.duration         = 3e-3;
    rfData.cycles           = 2;
    rfData.sliceThickness   = sliceThickness;
    if strcmpi(sliceSelectList{ss},'on')
        rfData.doSliceSelect = 1;
    else
        rfData.doSliceSelect = 0;
    end
    rfData.preRWScale    	= 0;
    rfData.postRWScale    	= 0;
    rfData.keepTimingSS     = 0;
    
    %% get the sequence
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
    
    
    %% run reference: analytical double at dt = 1e-7
    % get the sequence
    [ pulseSequenceRef ] = spinTwin.test.seq.singleRF( ...
        rfData, mrSystem, 1e-8);
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
    
    testNum = 0;
    for ee = 1:length(engineList)
        for pp = 1:length(precisionList)
            
            % case entries
            engine      = engineList{ee};
            odeMethod   = odeMethodList{ee};
            precision   = precisionList{pp};
            
            % update the fields
            simControl.simulationEngine = engine;
            simControl.odeMethod        = odeMethod;
            simControl.precision        = precision;
            
            %% call the kernel
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
                    [solSim, stats] = spinTwin.fwdBloch.runPhasor( ...
                        solSim, pulseSequence, simModel,...
                        motionModel, simControl, dbgControl );
                case lower('diffusion')
                    [solSim, stats] = spinTwin.fwdBloch.runDiffusion(...
                        solSim, pulseSequence, simModel,...
                        motionModel, simControl, dbgControl );
                otherwise % basic Bloch
                    [solSim, stats] = spinTwin.fwdBloch.runBloch(...
                        solSim, pulseSequence, simModel,...
                        motionModel, simControl, dbgControl );
            end
            
            fprintf('...');
            
            % magnitude and phase
            solSim.Mm = abs(solSim.Mx + 1j*solSim.My);
            solSim.Mp = angle(solSim.Mx + 1j*solSim.My);
            
            %% verify correctness with exponential
            errorMx = solSim.Mx-solRef.Mx;
            errorMy = solSim.My-solRef.My;
            errorMz = solSim.Mz-solRef.Mz;
            
            %% store test results
            testCase.engine     = engine;
            testCase.ode        = odeMethod;
            testCase.precision  = precision;
            testCase.errorMx    = errorMx;
            testCase.errorMy    = errorMy;
            testCase.errorMz    = errorMz;
            testCase.stats      = stats;
            % assign
            testNum = testNum + 1;
            testResults{testNum} = testCase;
            
            tKernel     = single(testCase.stats.tRF);
            numSteps    = single(testCase.stats.nRF);
            nVoxels     = single(testCase.stats.numVoxels);
            perf        = 1e12*tKernel/nVoxels/max(numSteps,1);
            
            %% plot results in different plots for each type
            figure(ss); hold on;
            perfAxes = {sprintf('%10s: \t %12s \t %10s',...
                testCase.engine, testCase.ode, testCase.precision)};
            barh(categorical(perfAxes),perf);
            xlabel('average ps/voxel/step');
            title(sprintf('RF Slice Selection %s, %d Voxels %d Time steps\n(%s, kernel %s)',...
                sliceSelectList{ss}, nVoxels, numSteps,...
                gpuDevice().Name, strrep(simControl.kernelVer,'_',' ')));
            grid on;
            
            
        end
    end
    
    fprintf('\n');
    fprintf('\n SLICE SELECTION %s, %d Voxels, %d time steps', sliceSelectList{ss}, stats.numVoxels, stats.nRF);
    fprintf('\n ENGINE    \t ODE      \t PRECISION   \t ERROR Mx \t ERROR My \t ERROR Mz \t PERFORMANCE');
    for ii = 1:testNum
        % get test
        testCase = testResults{ii};
        
        tKernel     = single(testCase.stats.tRF);
        numSteps    = single(testCase.stats.nRF);
        nVoxels     = single(testCase.stats.numVoxels);
        perf        = 1e12*tKernel/nVoxels/max(numSteps,1);
        
        % print results
        fprintf('\n %10s \t  %10s \t %10s \t %.3e \t %.3e \t %.3e \t %.1fps/voxel/step',...
            testCase.engine, testCase.ode, testCase.precision,...
            max(abs(testCase.errorMx)), max(abs(testCase.errorMy)), max(abs(testCase.errorMz)),...
            perf);
    end
    fprintf('\n');
    
end
