%
% SIMULATOR.EXAMPLE.RF180EXAMPLE
%
%	Script to verify the effect of a 180 refocusing pulse 
%   when initial phase (and its derivative) is not 0.
%   Generates a PGSE encoding (90 + GR + 180 + GR), and checks that 
%   the final magnetization phase and derivative is 0
%   The script runs analytical and phasor simulations, and compares the 
%   magnetizations at the end in the slice
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%
clear all; close all; clc;

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
doSliceSelect   = 0;

%% PGSE parameters
gradDir     = 'x';
TAU         = 50e-3;
gradAmp     = 0.030;
gradTime    = 20e-3;

%% initial magnetization state
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
t1Value = 10.500;
t2Value = 10.300;
switch lower(gradDir)
    case 'x'
        fovX = 0.100; dx = 5e-4;
    case 'y'
        fovY = 0.100; dy = 5e-4;
    case 'z'
        fovZ = 0.100; dz = 5e-4;
end

%% control for the experiment
expControl = data.experiment.initializeExpControl();
mrSystem = expControl.mrSystem;
motionModel = [];
expControl.debug.debugMode = 1;
expControl.sequence.dtRF   = 5e-6;
expControl.application = 'test';
expControl.simulation.kernelPtx = '/efs-mount-point/S20/PROJECTS/edutool/kernels/v26_20210128_ufm_sm37.ptx';

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


%% generate spin Model
[spinModel] = models.homogeneousCubeModel( fovX, fovY, fovZ, dx, dy, dz,...
    t1Value, t2Value, b0, mu, expControl.debug.debugMode);

% get sizes
numCoils = spinModel.slice{1}.model.numRxCoils;
numIso   = spinModel.slice{1}.model.numIsochromats;
numRxs   = pulseSequence.numRxs;


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

figure(1);
switch lower(gradDir)
    case 'x'
        position = spinModel.slice{1}.model.x;
    case 'y'
        position = spinModel.slice{1}.model.y;
    case 'z'
        position = spinModel.slice{1}.model.z;
end

figure(1);
subplot(5,1,1);  hold on;
plot(position, solAnalytical.Mx, 'LineWidth', 2 ); ylabel('Mx');
subplot(5,1,2);  hold on;
plot(position, solAnalytical.My, 'LineWidth', 2 ); ylabel('My');
subplot(5,1,3);  hold on;
plot(position, solAnalytical.Mz, 'LineWidth', 2 ); ylabel('Mz');
subplot(5,1,4);  hold on;
plot(position, abs(solAnalytical.Mm), 'LineWidth', 2 ); ylabel('|Mxy|');
subplot(5,1,5);  hold on;
plot(position, solAnalytical.Mp*180/pi, 'LineWidth', 2 ); ylabel('Phase(Mxy)');

figure(2); hold on;
plot(position, [diff(solAnalytical.Mp)/(position(2)-position(1));0], 'LineWidth', 2 ); ylabel('dPhase(Mxy)');


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

figure(1); hold on;        
subplot(5,1,1);  hold on;
plot(position, solPhasor.Mx, 'LineWidth', 2 ); 
ylabel('Mx'); legend('Analytic','Phasor');
subplot(5,1,2);  hold on;
plot(position, solPhasor.My, 'LineWidth', 2 ); 
ylabel('My'); legend('Analytic','Phasor');
subplot(5,1,3);  hold on;
plot(position, solPhasor.Mz, 'LineWidth', 2 ); 
ylabel('Mz'); legend('Analytic','Phasor');
subplot(5,1,4);  hold on;
plot(position, abs(solPhasor.Mm), 'LineWidth', 2 ); 
ylabel('|Mxy|'); legend('Analytic','Phasor');
subplot(5,1,5);  hold on;
plot(position, solPhasor.Mp*180/pi, 'LineWidth', 2 ); 
ylabel('Phase(Mxy)'); legend('Analytic','Phasor');
xlabel('time (s)'); 

switch lower(gradDir)
    case 'x'
        deriv = solPhasor.dMpDx;
    case 'y'
        deriv = solPhasor.dMpDy;
    case 'z'
        deriv = solPhasor.dMpDz;
end

figure(2); hold on;
plot(position, deriv, 'LineWidth', 2 ); 


%% call the Phasor kernel
expControl.simulation.simulationEngine = 'ODE';
expControl.simulation.kernelPtx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/v30_20210219_ufm_sm37.ptx';
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

figure(1); hold on;
subplot(5,1,1);  hold on;
plot(position, solPhasor.Mx, 'LineWidth', 2 ); 
ylabel('Mx'); legend('Analytic','Phasor');
subplot(5,1,2);  hold on;
plot(position, solPhasor.My, 'LineWidth', 2 ); 
ylabel('My'); legend('Analytic','Phasor');
subplot(5,1,3);  hold on;
plot(position, solPhasor.Mz, 'LineWidth', 2 ); 
ylabel('Mz'); legend('Analytic','Phasor');
subplot(5,1,4);  hold on;
plot(position, abs(solPhasor.Mm), 'LineWidth', 2 ); 
ylabel('|Mxy|'); legend('Analytic','Phasor');
subplot(5,1,5);  hold on;
plot(position, solPhasor.Mp*180/pi, 'LineWidth', 2 ); 
ylabel('Phase(Mxy)'); legend('Analytic','Phasor');
xlabel('time (s)'); 


switch lower(gradDir)
    case 'x'
        deriv = solPhasor.dMpDx;
    case 'y'
        deriv = solPhasor.dMpDy;
    case 'z'
        deriv = solPhasor.dMpDz;
end

figure(2); hold on;
plot(position, deriv, 'LineWidth', 2 ); 
ylabel('dPhase(Mxy)'); xlabel('position z (m)'); legend('Analytic','Phasor');

