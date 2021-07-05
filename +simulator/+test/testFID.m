clear all;
close all;
clc;

expControl.simulation.simulationEngine  = 'numerical'; % analytical / phasor / diffusion / numerical
expControl.simulation.precision         = 'single'; % single / double, only available in DF parts
expControl.simulation.odeMethod         = 'explicit'; % for numerical RF, explicit / adaptiveExp / implicit / adaptiveImp

%% define how to simulate
%   RF for uses specific kernel as defined for analytical/phasor/diffusion/numerical
%   GR for numerical/analytical uses same (analytical), for phasor/diffusion uses specific
%   RO for numerical/analytical uses same (analytical), for phasor/diffusion uses specific
%   DF for numerical/analytical uses same (analytical) with given precision, not available in phasor / diffusion
usePart = 'RF';

expControl.simulation.kernelPtx = '/efs-mount-point/S20/PROJECTS/edutool/kernels/v26_sm37.ptx'; %'/efs-mount-point/S20/PROJECTS/edutool/kernels/v26_sm37.ptx';
expControl.simulation.threads   = 256;
expControl.application          = 'test';
expControl.debug.debugMode      = 1;
expControl.debug.debugFile      = [];
expControl.debug.waitBarBE      = 0;
motionModel = [];
GPUindex    = 1;

%% define time and the values of T1 and T2
dt = 1e-6; 
totalTime = 0.2;
t1Val = 1e-1:1e-3:10;
t2Val = 0.5;

% create combination of values
T1 = kron(t1Val(:),ones(length(t2Val),1));
T2 = kron(ones(length(t1Val),1),t2Val(:));
numIso = length(T1);

%% define conditions and initial magnetization state
b0                  = 1.0;
mu                  = 1.0;
initialMagnitude    = b0*mu;
initialPhase        = 0.0;
initialFA           = 90.0;

%% create model
fovX = 0.010;
fovY = 1e-4;
fovZ = 1e-4;
dx   = 1e-4;
dy   = 1e-4;
dz   = 1e-4;

% make sure number of voxels is larger than T1/T2 values
nX = ceil(numIso/((fovY/dy)*(fovZ/dz)));
fovX = nX*dz;

% generate spin Model
[spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ,...
    dx, dy, dz, 1, 1, b0, mu);

%% modify the T1/T2 values to be different for each voxel
spinModel.slice{1}.model.numTissues             = numIso;
spinModel.slice{1}.model.numProperties          = 6;
spinModel.slice{1}.model.tissueValues           = zeros(numIso,6);
spinModel.slice{1}.model.tissueValues(:,1)      = T1(:);
spinModel.slice{1}.model.tissueValues(:,2)      = T2(:);
spinModel.slice{1}.model.tissueType(1:numIso)   = reshape(1:numIso,numIso,1);

% assign to simModel
simModel = spinModel.slice{1}.model;

%% generate a zeros pulse sequence
nSteps = ceil(totalTime/dt);
time = linspace(dt,totalTime,nSteps);

pulseSequence.name       = 'Dummy';
pulseSequence.type       = 'N/A';
pulseSequence.endEvent   = 'none'; % indicates what happens at the end
pulseSequence.totalTime  = time(end); % total time in seconds
pulseSequence.numSteps   = nSteps;   % number of time steps
pulseSequence.numRxs     = 1;   % number of readout points
pulseSequence.gamma      = 42.577478518e6; 

%% main sequence data
% waveforms
pulseSequence.time      = reshape(time,nSteps,1); % times
pulseSequence.rxSignal  = zeros(pulseSequence.numSteps,1); % receiver readout
pulseSequence.rxSignal(end) = 1;

pulseSequence.timeDiff  = reshape(diff([0; pulseSequence.time]),nSteps,1); % time deltas
pulseSequence.swcSignal = zeros(pulseSequence.numSteps,1); % software crusher
pulseSequence.gxSignal  = zeros(pulseSequence.numSteps,1); % x gradient
pulseSequence.gySignal  = zeros(pulseSequence.numSteps,1); % y gradient
pulseSequence.gzSignal  = zeros(pulseSequence.numSteps,1); % z gradient
pulseSequence.rfmSignal = zeros(pulseSequence.numSteps,1); % RF magnitude
pulseSequence.rfpSignal = zeros(pulseSequence.numSteps,1); % RF phase
pulseSequence.rffSignal = zeros(pulseSequence.numSteps,1); % RF frequency

% for splitting into parts (RF vs Non-RF)
pulseSequence.numParts     = 1; % number of parts
pulseSequence.partType{1}  = usePart; % type of part: RF / RO / GR / DF
pulseSequence.partLimits   = [1, nSteps]-1; % index start/end of part
% for indicating RO parts
pulseSequence.rxLimits     = [nSteps, nSteps];     

%% initialize solution
timeSolution = data.solution.initialize(...
        simModel.numIsochromats, pulseSequence.numRxs, simModel.numRxCoils);

%% initial magnetizations
% get sizes
numCoils = simModel.numRxCoils;
numIso   = simModel.numIsochromats;
numRxs   = pulseSequence.numRxs;

timeSolution.Mm         = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
timeSolution.Mp         = initialPhase*ones(numIso,1);
timeSolution.Mx         = timeSolution.Mm.*cos(initialPhase);
timeSolution.My         = timeSolution.Mm.*sin(initialPhase);
timeSolution.Mz         = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
timeSolution.dMpDx      = 0*timeSolution.Mp;
timeSolution.dMpDy      = 0*timeSolution.Mp;
timeSolution.dMpDz      = 0*timeSolution.Mp;
timeSolution.Sx         = zeros(numCoils*numRxs,1);
timeSolution.Sy         = zeros(numCoils*numRxs,1);
timeSolution.Sz         = zeros(numCoils*numRxs,1);
    
%% references: target result
E2 = exp(-pulseSequence.totalTime./simModel.tissueValues(simModel.tissueType(:),2));
E1 = exp(-pulseSequence.totalTime./simModel.tissueValues(simModel.tissueType(:),1));
refMm = timeSolution.Mm.*E2;
refMp = timeSolution.Mp;
refMz = initialMagnitude + ( timeSolution.Mz - initialMagnitude ).*E1;

%% simulate

switch lower(expControl.simulation.simulationEngine)
    case lower('phasor')
        [timeSolution] = simulator.bloch.kernel.runPhasor( ...
            timeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
    case lower('diffusion')
        [timeSolution] = simulator.diffusion.kernel.runSDW(...
            timeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
    case lower('analytical')
        [timeSolution] = simulator.bloch.kernel.runAnalytical(...
            timeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
    otherwise % 'numerical' ODE for RF, analytical for RO/GR and DF
        [timeSolution] = simulator.bloch.kernel.runNumerical(...
            timeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
end


%% plot
figure(1); plot(refMm(:)); hold on; 
plot(timeSolution.Mm(:));
legend('M_x_y reference', 'M_x_y simulated');
figure(2); plot(1000*simModel.tissueValues(:,1),refMz(:)); hold on; 
plot(1000*simModel.tissueValues(:,1),timeSolution.Mz(:));
legend('M_z reference', 'M_z simulated');
title(['Mz after ',num2str(pulseSequence.totalTime),'sec post-90RF for spin with given T1'])
xlabel('T1 value (ms)')
ylabel('Mz')

