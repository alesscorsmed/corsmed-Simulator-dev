%
% SIMULATOR.EXAMPLE.FLIPANGLEEXAMPLE
%
%	Script to verify the effect of a slice selective RF pulse.
%   The user can define the type of RF, flip angle, phase, slice selection,
%   pre and post rewind of slice selection gradient, discretization, ...
%   The script runs analytical and phasor simulations, and compares the 
%   magnetizations at the end in the slice
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

%% define RF conditions
gamma           = 42.577478518e6;
b0              = 1.0;
mu              = 1;
f0              = gamma*b0;
dt              = 1e-6;
sliceThickness  = 6e-3;
rfDuration      = 3e-3;
rfFlipAngle     = 90;
rfPhase         = 0.0;
rfCycles        = 2;
maxGStrenght    = 40e-3;
gradSlewrate    = 150.0;
preRWScale      = 0;
postRWScale     = 0;
doSliceSelect   = 1;

%% initial magnetization state
initialMagnitude    = b0*mu;
initialPhase        = 0.0;
initialFA           = 0.0;

%% domain
fovX = 0.001;
fovY = 0.001;  
fovZ = 2*sliceThickness;
dx   = 1e-3;
dy   = 1e-3;
dz   = 1e-4;
t1Value = 5.000;
t2Value = 3.000;

%% control for the experiment
[expControl] = data.expControl.initialize();
expControl.kernelPtx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_23.ptx';

%% mr System specs
[mrSystem] = data.mrSystem.initialize();
mrSystem.maxGStrenght   = maxGStrenght;
mrSystem.SlewRate       = gradSlewrate;

%% generate spin Model
[spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ, dx, dy, dz,...
    t1Value, t2Value, b0, mu);

%% generate the sequence
pulseSequence = data.pulseSequence.initialize();
%% assign values
pulseSequence.mainRF.flipAngle          = rfFlipAngle; % in degrees
pulseSequence.mainRF.phase              = rfPhase; % in rads
pulseSequence.mainRF.duration           = rfDuration;
pulseSequence.mainRF.cycles             = rfCycles;
pulseSequence.mainRF.sliceThickness     = sliceThickness;
pulseSequence.mainRF.doSliceSelect      = doSliceSelect;
pulseSequence.mainRF.preRWScale         = preRWScale;
pulseSequence.mainRF.postRWScale        = postRWScale;
%% populate waveforms
[pulseSequence]=sequence.dummy.singleRF(pulseSequence, mrSystem, expControl);

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

figure(1);
z = spinModel.slice{1}.model.z;
subplot(5,1,1);  hold on;
plot(z, solAnalytical.Mx, 'LineWidth', 2 ); ylabel('Mx');
subplot(5,1,2);  hold on;
plot(z, solAnalytical.My, 'LineWidth', 2 ); ylabel('My');
subplot(5,1,3);  hold on;
plot(z, solAnalytical.Mz, 'LineWidth', 2 ); ylabel('Mz');
subplot(5,1,4);  hold on;
plot(z, abs(solAnalytical.Mm), 'LineWidth', 2 ); ylabel('|Mxy|');
subplot(5,1,5);  hold on;
plot(z, solAnalytical.Mp*180/pi, 'LineWidth', 2 ); ylabel('Phase(Mxy)');


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

figure(1);
z = spinModel.slice{1}.model.z;
subplot(5,1,1);  hold on;
plot(z, solPhasor.Mx, 'LineWidth', 2 ); 
ylabel('Mx'); legend('Analytic','Phasor');
subplot(5,1,2);  hold on;
plot(z, solPhasor.My, 'LineWidth', 2 ); 
ylabel('My'); legend('Analytic','Phasor');
subplot(5,1,3);  hold on;
plot(z, solPhasor.Mz, 'LineWidth', 2 ); 
ylabel('Mz'); legend('Analytic','Phasor');
subplot(5,1,4);  hold on;
plot(z, abs(solPhasor.Mm), 'LineWidth', 2 ); 
ylabel('|Mxy|'); legend('Analytic','Phasor');
subplot(5,1,5);  hold on;
plot(z, solPhasor.Mp*180/pi, 'LineWidth', 2 ); 
ylabel('Phase(Mxy)'); legend('Analytic','Phasor');
xlabel('time (s)'); 

figure(2); clf; hold on;
plot(z, [diff(solAnalytical.Mp)/(z(2)-z(1));0], 'LineWidth', 2 ); ylabel('dPhase(Mxy)');
plot(z, solPhasor.dMpDz, 'LineWidth', 2 ); 
ylabel('dPhase(Mxy)'); xlabel('position z (m)'); legend('Analytic','Phasor');

figure(3); plot(z, solPhasor.Mm-M0z);