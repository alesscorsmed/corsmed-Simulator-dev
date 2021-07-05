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

clear all; close all; clc;

%% discretizations to test
numberVoxels    = [50, 25, 10, 5, 1];

%% initialize the sim Control with the default versioning for testing
[simControl] = spinTwin.setup.initializeSimControl();

%% prepare data
dt  = 1e-6;
% gradient and echo time
grData.echoTime         = 5e-3;
grData.duration      	= 2.5e-3;
grData.amplitude     	= 0.040;
grData.direction        = 'z'; % encoding direction 'x' / 'y' / 'z'
% mrSystem
mrSystem.b0             = 1.5000;
mrSystem.maxGStrenght   = 0.0400;
mrSystem.SlewRate       = 0.0;
% motion
motionModel             = [];
% dbg
dbgControl.mode         = 1;
dbgControl.file         = [];

%% initial magnetization state
b0                  = 1.0;
mu                  = 1.0;
initialMagnitude    = b0*mu;
initialPhase        = 0.0;
initialFA           = 90.0;

%% domain
fovX = 0.005;
fovY = 0.005;
fovZ = 0.005;
t1Value = 0.500;
t2Value = 0.250;

%% populate the sequence
[pulseSequence, grData] = spinTwin.test.seq.pgReversed( ...
    grData, mrSystem, dt);
                
figure(1); subplot(2,1,1); hold on;
plot(pulseSequence.time, pulseSequence.gzSignal, 'LineWidth',2 );
xlabel('time (s)'); ylabel('Gradient (T/m)');


%% loop on the pre-defined discretizations
for nVox = numberVoxels
    
    %% generate spin Model for given resolution
    dx = fovX/nVox;
    dy = fovY/nVox;
    dz = fovZ/nVox;
    [spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
        fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
    
    % shift the coordinates of the model for more impact of the gradient
    spinModel.slice{1}.model.r3D = spinModel.slice{1}.model.r3D + fovX;
    spinModel.slice{1}.model.x = spinModel.slice{1}.model.x + fovX;
    spinModel.slice{1}.model.y = spinModel.slice{1}.model.y + fovX;
    spinModel.slice{1}.model.z = spinModel.slice{1}.model.z + fovX;
    
    %% get dimensions and initialize magnetizations
    numCoils = spinModel.slice{1}.model.numRxCoils;
    numIso   = spinModel.slice{1}.model.numIsochromats;
    numRxs   = pulseSequence.numRxs;
    
    %% initial magnetizations
    M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
    M0p = initialPhase*ones(numIso,1);
    M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
    M0x = M0m.*cos(initialPhase);
    M0y = M0m.*sin(initialPhase);
    
    %% call the standard Bloch kernel, no signal integration
    solBloch = data.simulation.initializeSolution(numIso,numRxs,numCoils);
    solBloch.indexes    = [];
    solBloch.Mx         = 1*M0x;
    solBloch.My         = 1*M0y;
    solBloch.Mz         = 1*M0z;
    solBloch.dMpDx      = 0*M0p;
    solBloch.dMpDy      = 0*M0p;
    solBloch.dMpDz      = 0*M0p;
    solBloch.Sx         = zeros(numCoils*numRxs,1);
    solBloch.Sy         = zeros(numCoils*numRxs,1);
    solBloch.Sz         = zeros(numCoils*numRxs,1);
    % standard w/ no signal integration
    simControl.simulationEngine = 'Bloch';
    % kernel call
    [solBloch] = spinTwin.fwdBloch.runBloch(...
        solBloch, pulseSequence,spinModel.slice{1}.model,...
        motionModel, simControl, dbgControl );
    % magnitude and phase
    solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
    solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);
    
    figure(1); subplot(2,1,2); hold on;
    plot(pulseSequence.time, solBloch.Sx, 'LineWidth',2 );
    
    if nVox == 1
        
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
        % standard w/ no signal integration
        simControl.simulationEngine = 'Phasor';
        % kernel call
        [solPhasor] = spinTwin.fwdBloch.runPhasor(...
            solPhasor, pulseSequence,spinModel.slice{1}.model,...
            motionModel, simControl, dbgControl );
        % select correct derivative
        switch lower(grData.direction)
            case 'x'
                deriv = solPhasor.dMpDx;
            case 'y'
                deriv = solPhasor.dMpDy;
            case 'z'
                deriv = solPhasor.dMpDz;
        end
        % magnitude and phase
        solPhasor.Mm = abs(solPhasor.Mx + 1j*solPhasor.My);
        solPhasor.Mp = angle(solPhasor.Mx + 1j*solPhasor.My);
        
        figure(1); subplot(2,1,2); hold on;
        plot(pulseSequence.time, solPhasor.Sx, 'k--', 'LineWidth',2 );
        
    end
    
end

legend('50x50x50, 0.1mm', '25x25x25, 0.2mm', '10x10x10, 0.5mm', '5x5x5, 1.0mm', '1x1x1, 5.0mm', 'Phasor 1x1x1, 5.0mm');
xlabel('time (s)'); ylabel('Sx');
