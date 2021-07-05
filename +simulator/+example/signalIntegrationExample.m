%
% SIMULATOR.EXAMPLE.SIGNALINTEGRATIONEXAMPLE
%
%	Script to verify the effect of signal integration.
%   A Forward and Rewind gradient is applied, with some time in between.
%   A 5mm cube is discretized with different resolutions, and the 
%   analytical simulation is compated with the phasor simulation.
%   Phasor simulation produces the correct signal (correct integration)
%   for 5mm resolution (single voxel), whereas the analytical generates 
%   spurious echoes for coarse discretizations.
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

%% discretizations to test
numberVoxels    = [50, 25, 10, 5, 1];

%% define RF conditions
gamma           = 42.577478518e6;
b0              = 1.0;
mu              = 1;
f0              = gamma*b0;
dt              = 1e-6;

%% initial magnetization state
initialMagnitude    = b0*mu;
initialPhase        = 0.0;
initialFA           = 90.0;

%% domain
fovX = 0.005;
fovY = 0.005;
fovZ = 0.005;
t1Value = 0.500;
t2Value = 0.250;

%% control for the experiment
[expControl] = data.expControl.initialize();
expControl.kernelPtx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_23.ptx';

%% mr System specs
[mrSystem] = data.mrSystem.initialize();
mrSystem.maxGStrenght   = maxGStrenght;
mrSystem.SlewRate       = gradSlewrate;

%% generate the sequence
TAU = 5e-3; gradAmp = 0.03; gradTime = 2e-3; gradDir = 'z';
%% assign values
pulseSequence = data.pulseSequence.initialize();
pulseSequence.encPG.Dir = gradDir;
pulseSequence.encPG.TG  = gradTime;
pulseSequence.encPG.AG  = gradAmp;
pulseSequence.encPG.TAU = TAU;
%% populate the sequence
[pulseSequence] = sequence.pgse.reversedPG(...
    pulseSequence, mrSystem, expControl);
                
figure(1); subplot(2,1,1); hold on;
plot(pulseSequence.time, pulseSequence.gzSignal, 'LineWidth',2 );
xlabel('time (s)'); ylabel('Gradient (T/m)');


%% loop on the pre-defined discretizations
for nVox = numberVoxels
    
    %% generate spin Model for given resolution
    dx = fovX/nVox;
    dy = fovY/nVox;
    dz = fovZ/nVox;
    [spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ, dx, dy, dz,...
        t1Value, t2Value, b0, mu);
    
    % shift the coordinates of the model for more impact of the gradient
    spinModel.slice{1}.model.r3D = spinModel.slice{1}.model.r3D + fovX;
    spinModel.slice{1}.model.x = spinModel.slice{1}.model.x + fovX;
    spinModel.slice{1}.model.y = spinModel.slice{1}.model.y + fovX;
    spinModel.slice{1}.model.z = spinModel.slice{1}.model.z + fovX;
    
    %% initial magnetizations
    % get sizes
    numCoils = spinModel.slice{1}.model.numCoils;
    numIso   = spinModel.slice{1}.model.numIsochromats;
    numRxs   = pulseSequence.numRxs;
    
    M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
    M0p = initialPhase*ones(numIso,1);
    M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
    M0x = M0m.*cos(initialPhase);
    M0y = M0m.*sin(initialPhase);
    
    %% call the standard analytical kernel
    solAnalytical = data.solution.initialize(numIso,numRxs,numCoils);
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
    [solAnalytical] = simulator.bloch.kernel.runAnalytical(solAnalytical,...
        pulseSequence,spinModel.slice{1}.model,expControl);
    
    figure(1); subplot(2,1,2); hold on;
    plot(pulseSequence.time, solAnalytical.Sx, 'LineWidth',2 );
    
    if nVox == 1
        
        %% call the Phasor kernel
        solPhasor = data.solution.initialize(numIso,numRxs,numCoils);
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
        [solPhasor] = simulator.bloch.kernel.runPhasor(solPhasor,...
            pulseSequence,spinModel.slice{1}.model,expControl);
        
        figure(1); subplot(2,1,2); hold on;
        plot(pulseSequence.time, solPhasor.Sx, 'k--', 'LineWidth',2 );
        
    end
    
end

legend('50x50x50, 0.1mm', '25x25x25, 0.2mm', '10x10x10, 0.5mm', '5x5x5, 1.0mm', '1x1x1, 5.0mm', 'Phasor 1x1x1, 5.0mm');
xlabel('time (s)'); ylabel('Sx');
