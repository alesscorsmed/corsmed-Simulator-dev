%
% SPINTWIN.TEST.PERFORMANCE.FID
%
%	Runs an FID, compares results and performance
%
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

% cases to run
precisionList        = {'single', 'double'};
odeMethodList        = {'analytical', 'explicit'};
engineList           = {'bloch', 'phasor', 'diffusion'};
partTypeList         = {'RF', 'GR', 'RO'};
gzAmplitudeList      = {0.010};

%% general info for the test
% 10 ms relaxation at 1us discretization
dt          = 1e-6;
numSteps    = 1e5;
t2Limits    = [0.010,  5.000];
t1Limits    = [0.040,  5.000];
numT2       = 256;
numT1       = 64;
%% initial magnetization state
b0          = 1.0;
mu          = 1.0;
initialM    = b0*mu;
initialP    = 0.0;
initialFA   = 90.0; % in degrees

%% prepare the data required
% mrSystem
mrSystem.b0             = 1.5000;
mrSystem.maxGStrenght   = 0.0400;
mrSystem.SlewRate       = 0.0;
% motion
motionModel             = [];
% dbg
dbgControl.mode         = 1;
dbgControl.file         = [];

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
pulseSequence.gzSignal(:)   = 0.0;
pulseSequence.rxSignal(:)   = 1;
pulseSequence.totalTime     = dt*numSteps; % total time in seconds
pulseSequence.numSteps      = numSteps;   % number of time steps
pulseSequence.numRxs        = numSteps;   % number of readout points
% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts      = 1; % number of parts
pulseSequence.partType{1}   = 'RF'; % type of part: RF / RO / GR / DF
pulseSequence.partLimits    = [1, numSteps]; % index start/end of part
% for indicating RO parts
pulseSequence.rxLimits      = [numSteps, numSteps];

%% model
dx = 1e-3; dy = 1e-3; dz = 1e-3;
fovX = numT2*dx;
fovY = numT1*dy;
fovZ = 11*dz;
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

testNum = 0;
for gzAmplitude = cell2mat(gzAmplitudeList)

    %% update pulse sequence 
    pulseSequence.gzSignal(:)   = gzAmplitude;

    %% references: target result
    % cumulated area
    theta = sum(pulseSequence.timeDiff(:).*pulseSequence.gzSignal(:));
    R  = -2*pi*gamma*theta*simModel.z; % phase for each voxel
    E2 = exp(-pulseSequence.totalTime./simModel.tissueValues(simModel.tissueType(:),2));
    E1 = exp(-pulseSequence.totalTime./simModel.tissueValues(simModel.tissueType(:),1));
    M0m   = abs(M0x + 1j*M0y);
    refMx = M0m.*E2.*cos(R);
    refMy = M0m.*E2.*sin(R);
    refMm = abs(refMx + 1j*refMy);
    refMp = angle(refMx + 1j*refMy);
    refMz = initialM + ( M0z - initialM ).*E1;

    for tt = 1:length(partTypeList)
        for ee = 1:length(engineList)
            for oo = 1:length(odeMethodList)
                for pp = 1:length(precisionList)
                    
                    % case entries
                    partType    = partTypeList{tt};
                    engine      = engineList{ee};
                    odeMethod   = odeMethodList{oo};
                    precision   = precisionList{pp};
                    
                    
                    % update the fields
                    pulseSequence.partType{1}   = partType; % type of part: RF / RO / GR / DF
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
                    
                    % magnitude and phase
                    solSim.Mm = abs(solSim.Mx + 1j*solSim.My);
                    solSim.Mp = angle(solSim.Mx + 1j*solSim.My);
                    
                    %% verify correctness with exponential
                    errorMx = solSim.Mx-refMx;
                    errorMy = solSim.My-refMy;
                    errorMz = solSim.Mz-refMz;
                                        
                    %% prepare data to store and plot (if useful)
                    switch partType
                        case 'RF'
                            tKernel     = single(stats.tRF);
                            numSteps    = single(stats.nRF);
                            storeData   = 1;
                        case 'RO'
                            tKernel     = single(stats.tRO);
                            numSteps    = single(stats.nRO);
                            if strcmpi(odeMethod, 'analytical')
                                storeData = 1;
                            else % RO explicit defaults to analytical
                                storeData = 0; 
                            end
                        case 'GR'
                            tKernel     = single(stats.tGR);
                            numSteps    = single(stats.nGR);
                            if strcmpi(odeMethod, 'analytical')
                                storeData = 1;
                            else % GR explicit defaults to analytical
                                storeData = 0;
                            end
                        otherwise
                            tKernel     = single(stats.tDF);
                            numSteps    = single(stats.nDF);
                            if strcmpi(odeMethod, 'analytical')
                                storeData = 1;
                            else % DF explicit defaults to analytical
                                storeData = 0;
                            end
                    end
                    nVoxels = single(stats.numVoxels);
                    perf    = 1e12*tKernel/nVoxels/max(numSteps,1);
                    
                    %% store test results
                    if storeData
                        testCase.engine     = engine;
                        testCase.ode        = odeMethod;
                        testCase.part       = partType;
                        testCase.precision  = precision;
                        testCase.errorMx    = errorMx;
                        testCase.errorMy    = errorMy;
                        testCase.errorMz    = errorMz;
                        testCase.stats      = stats;
                        % assign
                        testNum = testNum + 1;
                        testResults{testNum} = testCase;
                        
                        %% plot results in different plots for each type
                        figure(tt); hold on;
                        perfAxes = {sprintf('%10s: \t %12s \t %10s',...
                            testCase.engine, testCase.ode, testCase.precision)};
                        barh(categorical(perfAxes),perf);
                        xlabel('average ps/voxel/step');
                        title(sprintf('FID case: %s type, %d Voxels %d Time steps\n(%s, kernel %s)',...
                            partType, nVoxels, numSteps,...
                            gpuDevice().Name, strrep(simControl.kernelVer,'_',' ')));
                        grid on;
                    end

                    
                    
                end
            end
        end
         
    end
end
                
fprintf('\n');
fprintf('\n ENGINE    \t ODE      \t PRECISION \t PART    \t ERROR Mx \t ERROR My \t ERROR Mz \t PERFORMANCE');
for ii = 1:testNum
    % get test
    testCase = testResults{ii};
    switch testCase.part
        case 'RF'
            tKernel     = single(testCase.stats.tRF);
            numSteps    = single(testCase.stats.nRF);
        case 'RO'
            tKernel     = single(testCase.stats.tRO);
            numSteps    = single(testCase.stats.nRO);
        case 'GR'
            tKernel     = single(testCase.stats.tGR);
            numSteps    = single(testCase.stats.nGR);
        otherwise
            tKernel     = single(testCase.stats.tDF);
            numSteps    = single(testCase.stats.nDF);
    end
    nVoxels = single(testCase.stats.numVoxels);
    perf    = 1e12*tKernel/nVoxels/max(numSteps,1);

    % print results
    fprintf('\n %10s \t  %10s \t %10s \t %10s \t %.3e \t %.3e \t %.3e \t %.1fps/voxel/step',...
       testCase.engine, testCase.ode, testCase.precision, testCase.part,...
       max(abs(testCase.errorMx)), max(abs(testCase.errorMy)), max(abs(testCase.errorMz)),...
       perf);   
end
fprintf('\n');



