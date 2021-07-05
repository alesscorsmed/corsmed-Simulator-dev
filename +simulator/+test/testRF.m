%clear all;
%close all;
clc;

expControl.simulation.simulationEngine  = 'numerical'; % analytical / phasor / diffusion / numerical
expControl.simulation.precision         = 'double'; % single / double, only available in DF parts
expControl.simulation.odeMethod         = 'adaptiveExp'; % for numerical RF, explicit / adaptiveExp / implicit / adaptiveImp

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
expControl.sequence.dtRF = 1e-6;
sliceThickness = 5e-3;
fliAngle = 180;
t1Val = ones(1024,1); %1e-1:1e-3:10;
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
initialFA           = 0.0;

%% create model
fovX = 1e-4;
fovY = 1e-4;
dx   = 1e-4;
dy   = 1e-4;

% make sure number of voxels is larger than T1/T2 values
nZ      = ceil(numIso/((fovX/dx)*(fovY/dy)));
fovZ    = 2*sliceThickness;
dz      = fovZ/nZ;

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

%% generate a RF slice selection sequence

%% mr System specs
[mrSystem] = data.mrSystem.initialize();
mrSystem.maxGStrenght   = 40e-3;
mrSystem.SlewRate       = 0.0;

%% generate the sequence
% rf info
mainRF                  = data.acquisition.initializeRF(fliAngle);
rfPulse.phase           = -pi/2; % in radians
mainRF.sliceThickness   = sliceThickness;
% assign values to avoid pre and post SS gradient
mainRF.doSliceSelect    = 1;
mainRF.preRWScale       = 0;
mainRF.postRWScale      = 0;
mainRF.keepTimingSS     = 0;

%% populate waveforms
[pulseSequence]=sequence.dummy.singleRF(mainRF, mrSystem, expControl);
pulseSequence.partType{1} = usePart;
pulseSequence.partLimits  = pulseSequence.partLimits - 1;

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
figure(1);
subplot(3,1,1); hold on
plot(simModel.z(:), timeSolution.Mx(:));
legend('M_x simulated');
xlabel('z (m)')
ylabel('M_x')
subplot(3,1,2); hold on
plot(simModel.z(:), timeSolution.My(:));
legend('M_y simulated');
xlabel('z (m)')
ylabel('M_y')
subplot(3,1,3); hold on
plot(simModel.z(:), timeSolution.Mz(:));
legend('M_z simulated');
xlabel('z (m)')
ylabel('M_z')

