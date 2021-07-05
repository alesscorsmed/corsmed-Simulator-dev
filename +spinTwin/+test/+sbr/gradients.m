%
% SPINTWIN.TEST.EXAMPLE.GRADIENTS
%
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

%% initialize the sim Control with the default versioning for testing
[simControl] = spinTwin.setup.initializeSimControl();

precision = 'double';
seqType   = 'fid';

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


%% model and data
b0          = 1.0;
mu          = 1.0;
t2Limits    = [0.010,  5.000];
t1Limits    = [0.040,  5.000];
numT2       = 256;
numT1       = 64;

%% model
dx = 1e-3; dy = 1e-3; dz = 1e-3;
fovX = numT2*dx;
fovY = numT1*dy;
fovZ = 3*dz;
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

% get dimensions
numCoils = simModel.numRxCoils;
numIso   = simModel.numIsochromats;

%% generate SBR model
sbrModel.numIsochromats = simModel.numIsochromats;
sbrModel.b0          	= simModel.b0;
sbrModel.mu             = simModel.mu;
sbrModel.x              = simModel.x;
sbrModel.y              = simModel.y;
sbrModel.z              = simModel.z;
sbrModel.resolution     = simModel.resolution;
sbrModel.rxCoilMapsX    = simModel.rxCoilMapsX;
sbrModel.rxCoilMapsY    = simModel.rxCoilMapsY;
sbrModel.numRxCoils     = simModel.numRxCoils;
sbrModel.pr             =  ones(numIso,1);
sbrModel.pi             = zeros(numIso,1);
sbrModel.r1             = zeros(numIso,1);
sbrModel.r2             = zeros(numIso,1);
sbrModel.bi             = zeros(numIso,1);
sbrModel.cs             = zeros(numIso,1);
% assign values
sbrModel.pr(:) = simModel.pd(:);
sbrModel.r1(:) = 1./t1Mat(:);   sbrModel.r1(isnan(sbrModel.r1)) = 1e4;
sbrModel.r2(:) = 1./t2Mat(:);   sbrModel.r2(isnan(sbrModel.r1)) = 1e5;

%% initial magnetizations
initialMagnitude = b0*mu;
initialPhase 	 = 0.0;
initialFA      	 = 90.0;
M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
M0p = initialPhase*ones(numIso,1);
M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
M0x = M0m.*cos(initialPhase);
M0y = M0m.*sin(initialPhase);


%% prepare sequence data
if strcmpi(seqType, 'pgse')
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
    %% get the sequence
    [ pulseSequence, grData ] = spinTwin.test.seq.pgse( ...
        grData, rfData, mrSystem, dt);
else
    %% get the FID sequence
    dt  = 1e-6;
    totalTime = 1e-2;
    [ pulseSequence ] = spinTwin.test.seq.fid(totalTime,dt);
end
numRxs   = pulseSequence.numRxs;


%% call the forward kernel
simControl.simulationEngine = 'Phasor';
simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
simControl.precision        = precision;
% initialize solution
solFWD = data.simulation.initializeSolution(numIso,numRxs,numCoils);
solFWD.indexes    = [];
solFWD.Mx         = 1*M0x;
solFWD.My         = 1*M0y;
solFWD.Mz         = 1*M0z;
solFWD.dMpDx      = 0*M0x;
solFWD.dMpDy      = 0*M0x;
solFWD.dMpDz      = 0*M0x;
solFWD.Sx         = zeros(numCoils*numRxs,1);
solFWD.Sy         = zeros(numCoils*numRxs,1);
solFWD.Sz         = zeros(numCoils*numRxs,1);
% kernel call
[solFWD, statsFWD] = spinTwin.fwdBloch.runPhasor(...
    solFWD, pulseSequence, simModel,...
    motionModel, simControl, dbgControl );

% plot
figure(1); hold on;        
subplot(3,1,1);  hold on;
plot(solFWD.Mx, 'LineWidth', 2 ); 
ylabel('Mx');
subplot(3,1,2);  hold on;
plot(solFWD.My, 'LineWidth', 2 ); 
ylabel('My');
subplot(3,1,3);  hold on;
plot(solFWD.Mz, 'LineWidth', 2 ); 
ylabel('Mz');


%% call the SBR kernel
pulseSequence.partType{1}   = 'RF'; % check deriv for type of part: RF / RO / GR
simControl.simulationEngine = 'Phasor';
simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
simControl.precision        = precision;
% initialize solution
solSBR = data.simulation.initializeSolution(numIso,numRxs,numCoils);
solSBR.indexes    = [];
solSBR.Mx         = 1*M0x;
solSBR.My         = 1*M0y;
solSBR.Mz         = 1*M0z;
solSBR.dPx        = 0*M0x;
solSBR.dPy        = 0*M0x;
solSBR.dPz        = 0*M0x;
solSBR.dR1x       = 0*M0x;
solSBR.dR1y       = 0*M0x;
solSBR.dR1z       = 0*M0x;
solSBR.dR2x       = 0*M0x;
solSBR.dR2y       = 0*M0x;
solSBR.dR2z       = 0*M0x;
solSBR.Sx         = zeros(numCoils*numRxs,1);
solSBR.Sy         = zeros(numCoils*numRxs,1);
solSBR.Srefx      = solFWD.Sx;
solSBR.Srefy      = solFWD.Sy;
% kernel call
[solSBR, statsSBR] = spinTwin.sbrBloch.runGradients(...
    solSBR, pulseSequence, sbrModel,...
    simControl, dbgControl );


% plot
figure(1); hold on;        
subplot(3,1,1);  hold on;
plot(solSBR.Mx, 'LineWidth', 2 ); 
ylabel('Mx');
subplot(3,1,2);  hold on;
plot(solSBR.My, 'LineWidth', 2 ); 
ylabel('My');
subplot(3,1,3);  hold on;
plot(solSBR.Mz, 'LineWidth', 2 ); 
ylabel('Mz');

figure(2); hold on;        
subplot(3,1,1);  hold on;
plot(-totalTime*solFWD.Mx, 'LineWidth', 2 ); 
plot(solSBR.dR2x, 'LineWidth', 2 ); 
ylabel('dR2x'); legend('Ref', 'Sim');
subplot(3,1,2);  hold on;
plot(-totalTime*solFWD.My, 'LineWidth', 2 ); 
plot(solSBR.dR2y, 'LineWidth', 2 ); 
ylabel('dR2y');  legend('Ref', 'Sim');
subplot(3,1,3);  hold on;
plot(solSBR.dR2z, 'LineWidth', 2 ); 
ylabel('dR2z');

figure(3); hold on;        
subplot(3,1,1);  hold on;
plot(solSBR.dR1x, 'LineWidth', 2 ); 
ylabel('dR1x');
subplot(3,1,2);  hold on;
plot(solSBR.dR1y, 'LineWidth', 2 ); 
ylabel('dR1y');
subplot(3,1,3);  hold on;
plot(-totalTime*(solFWD.Mz - initialMagnitude), 'LineWidth', 2 ); 
plot(solSBR.dR1z, 'LineWidth', 2 );
ylabel('dR1z');  legend('Ref', 'Sim');

figure()
ax1 = subplot(2,1,1);  hold on;
plot(-totalTime*solFWD.Mx, 'LineWidth', 2 ); 
plot(solSBR.dR2x, '--', 'LineWidth', 2 ); 
ylabel('dR2x'); legend('Ref', 'Sim');
ax2 = subplot(2,1,2);  hold on;
plot(-totalTime*(solFWD.Mz - initialMagnitude), 'LineWidth', 2 ); 
plot(solSBR.dR1z, '--','LineWidth', 2 );
ylabel('dR1z');  legend('Ref', 'Sim');
linkaxes([ax1, ax2], 'x');
