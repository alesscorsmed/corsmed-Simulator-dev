%
% SPINTWIN.TEST.EXAMPLE.PHASEREFOCUS
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

%% initialize the sim Control with the default versioning for testing
[simControl] = spinTwin.setup.initializeSimControl();

%% prepare data
dt  = 1e-6;
% gradient and echo time
grData.echoTime         = 50e-3;
grData.duration      	= 25e-3;
grData.amplitude     	= 0.040;
grData.direction        = 'x'; % encoding direction 'x' / 'y' / 'z'
% RF
rfData.type             = 'sinc';
rfData.flipAngle      	= 90; % in degrees
rfData.phase         	= 0.0; % in rads
rfData.duration         = 1e-3;
rfData.cycles           = 2;
rfData.sliceThickness   = 3e-3;
rfData.doSliceSelect 	= 0;
rfData.preRWScale    	= 0;
rfData.postRWScale    	= 0;
rfData.keepTimingSS     = 0;
% mrSystem
mrSystem.b0             = 1.5000;
mrSystem.maxGStrenght   = 0.0400;
mrSystem.SlewRate       = 0.0;
% motion
motionModel             = [];
% dbg
dbgControl.mode         = 1;
dbgControl.file         = [];

%% get the sequence
[ pulseSequence, grData ] = spinTwin.test.seq.pgse( ...
    grData, rfData, mrSystem, dt);

%% domain and model
b0 	 = 1.0;
mu 	 = 1.0;
fovX = 0.001;
fovY = 0.001;  
fovZ = 0.001;
dx   = 1e-3;
dy   = 1e-3;
dz   = 1e-3;
t1Value = 10.500;
t2Value = 10.300;
% discretize on the gradient direction
switch lower(grData.direction)
    case 'x'
        fovX = 0.100; dx = 5e-4;
    case 'y'
        fovY = 0.100; dy = 5e-4;
    case 'z'
        fovZ = 0.100; dz = 5e-4;
end
% generate spin Model
[spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
    fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );
% voxel positions
switch lower(grData.direction)
    case 'x'
        position = spinModel.slice{1}.model.x;
    case 'y'
        position = spinModel.slice{1}.model.y;
    case 'z'
        position = spinModel.slice{1}.model.z;
end

%% get dimensions and initialize magnetizations
numCoils = spinModel.slice{1}.model.numRxCoils;
numIso   = spinModel.slice{1}.model.numIsochromats;
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

%% call the standard analytical kernel
simControl.simulationEngine = 'Bloch'; 
simControl.odeMethod        = 'analytical'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
% initialize solution
solBloch = data.simulation.initializeSolution(numIso,numRxs,numCoils);
solBloch.indexes    = [];
solBloch.Mx         = 1*M0x;
solBloch.My         = 1*M0y;
solBloch.Mz         = 1*M0z;
solBloch.dMpDx      = 0*M0x;
solBloch.dMpDy      = 0*M0x;
solBloch.dMpDz      = 0*M0x;
solBloch.Sx         = zeros(numCoils*numRxs,1);
solBloch.Sy         = zeros(numCoils*numRxs,1);
solBloch.Sz         = zeros(numCoils*numRxs,1);
% kernel call
[solBloch] = spinTwin.fwdBloch.runBloch(...
    solBloch, pulseSequence,spinModel.slice{1}.model,...
    motionModel, simControl, dbgControl );
% magnitude and phase
solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);

% plot
figure(1);
subplot(5,1,1);  hold on;
plot(position, solBloch.Mx, 'LineWidth', 2 ); ylabel('Mx');
subplot(5,1,2);  hold on;
plot(position, solBloch.My, 'LineWidth', 2 ); ylabel('My');
subplot(5,1,3);  hold on;
plot(position, solBloch.Mz, 'LineWidth', 2 ); ylabel('Mz');
subplot(5,1,4);  hold on;
plot(position, abs(solBloch.Mm), 'LineWidth', 2 ); ylabel('|Mxy|');
subplot(5,1,5);  hold on;
plot(position, solBloch.Mp*180/pi, 'LineWidth', 2 ); ylabel('Phase(Mxy)');
figure(2); hold on;
plot(position, [diff(solBloch.Mp)/(position(2)-position(1));0], 'LineWidth', 2 ); ylabel('dPhase(Mxy)');


%% call the Phasor kernel (Analytical)
simControl.simulationEngine = 'Phasor'; 
simControl.odeMethod        = 'analytical'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
% initialize solution
solPhasor = data.simulation.initializeSolution(numIso,numRxs,numCoils);
solPhasor.indexes    = [];
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
% plot
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
figure(2); hold on;
plot(position, deriv, 'LineWidth', 2 ); 


%% call the Phasor kernel (explicit ODE)
simControl.simulationEngine = 'Phasor'; 
simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
solPhasorODE = data.simulation.initializeSolution(numIso,numRxs,numCoils);
solPhasorODE.indexes    = [];
solPhasorODE.Mm         = 1*M0m;
solPhasorODE.Mp         = 1*M0p;
solPhasorODE.Mx         = 1*M0x;
solPhasorODE.My         = 1*M0y;
solPhasorODE.Mz         = 1*M0z;
solPhasorODE.dMpDx      = 0*M0p;
solPhasorODE.dMpDy      = 0*M0p;
solPhasorODE.dMpDz      = 0*M0p;
solPhasorODE.Sx         = zeros(numCoils*numRxs,1);
solPhasorODE.Sy         = zeros(numCoils*numRxs,1);
solPhasorODE.Sz         = zeros(numCoils*numRxs,1);
% kernel call
[solPhasorODE] = spinTwin.fwdBloch.runPhasor(...
    solPhasorODE, pulseSequence,spinModel.slice{1}.model,...
    motionModel, simControl, dbgControl );
% select correct derivative
switch lower(grData.direction)
    case 'x'
        deriv = solPhasorODE.dMpDx;
    case 'y'
        deriv = solPhasorODE.dMpDy;
    case 'z'
        deriv = solPhasorODE.dMpDz;
end
% magnitude and phase
solPhasorODE.Mm = abs(solPhasorODE.Mx + 1j*solPhasorODE.My);
solPhasorODE.Mp = angle(solPhasorODE.Mx + 1j*solPhasorODE.My);

figure(1); hold on;
subplot(5,1,1);  hold on;
plot(position, solPhasorODE.Mx, 'LineWidth', 2 ); 
ylabel('Mx'); 
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
subplot(5,1,2);  hold on;
plot(position, solPhasorODE.My, 'LineWidth', 2 ); 
ylabel('My'); 
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
subplot(5,1,3);  hold on;
plot(position, solPhasorODE.Mz, 'LineWidth', 2 ); 
ylabel('Mz'); 
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
subplot(5,1,4);  hold on;
plot(position, abs(solPhasorODE.Mm), 'LineWidth', 2 ); 
ylabel('|Mxy|'); 
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
subplot(5,1,5);  hold on;
plot(position, solPhasorODE.Mp*180/pi, 'LineWidth', 2 ); 
ylabel('Phase(Mxy)'); 
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
xlabel('time (s)'); 
figure(2); hold on;
plot(position, deriv, 'LineWidth', 2 ); 
ylabel('dPhase(Mxy)'); xlabel('position z (m)'); 
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');

