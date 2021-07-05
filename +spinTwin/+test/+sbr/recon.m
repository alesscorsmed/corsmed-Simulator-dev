%
% SPINTWIN.SBR.EXAMPLE.RECON
%
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

%% prepare the experiment
% sizes for the original acquisition
dx = 1e-3;      dy = 1e-3;      dz = 1e-3;
fovX = 0.100;   fovY = 0.100;   fovZ = 3*dz;
% generate target model and sequence
[pulseSequence,simModel,motionModel,mrSystem] = ...
    spinTwin.test.sbr.utils.prepareExperiment( fovX, fovY, fovZ, dx, dy, dz );

%% get dimensions
numCoils = simModel.numRxCoils;
numIso   = simModel.numIsochromats;
numRxs   = pulseSequence.numRxs;
%% initial magnetizations
initialMagnitude = simModel.b0*simModel.mu;
initialPhase 	 = 0.0;
initialFA      	 = 0.0;
M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
M0p = initialPhase*ones(numIso,1);
M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
M0x = M0m.*cos(initialPhase);
M0y = M0m.*sin(initialPhase);


%% initialize the sim Control with the default versioning for testing
[simControl] = spinTwin.setup.initializeSimControl();
% dbg
dbgControl.mode = 1;
dbgControl.file = [];

%% call the forward kernel to generate the 'acquired' data
simControl.simulationEngine = 'Phasor';
simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
simControl.precision        = 'double';
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

% we can do standard FFT
kSpace = transpose(reshape(solFWD.Sx + 1j*solFWD.Sy, numCoils, numRxs));
kSpace = reshape(kSpace,[],pulseSequence.numEnc);
iSpace = fftshift(ifftn(ifftshift(kSpace)));


%% prepare the data required for the SBR
initialPD = 1.0;
initialT1 = 1e3;
initialT2 = 1e4;
reconFOVx = fovX;
reconFOVy = fovY;
reconFOVz = fovZ;
reconNx   = 32;
reconNy   = 32;
reconNz   = 1;

%% generate base model
[sbrModel] = spinTwin.test.sbr.utils.generateSBRmodel( ...
    reconFOVx, reconFOVy, reconFOVz, reconNx, reconNy, reconNz,...
    initialPD, initialT1, initialT2, simModel.b0, simModel.mu);

%% generate base solution
numReconIso = sbrModel.numIsochromats;
% initialize solution
solSBR.indexes = [];
% initial conditions
solSBR.initialMagnitude = simModel.b0*simModel.mu;
solSBR.initialPhase 	= 0.0;
solSBR.initialFA      	= 0.0;
M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numReconIso,1);
M0p = initialPhase*ones(numReconIso,1);
M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numReconIso,1);
M0x = M0m.*cos(initialPhase);
M0y = M0m.*sin(initialPhase);
% Magnetizations at the center of the voxel
solSBR.Mx = M0x;
solSBR.My = M0y;
solSBR.Mz = M0z;
% Spatial derivatives of the Phase
solSBR.dPx = zeros(numReconIso,1);
solSBR.dPy = zeros(numReconIso,1);
solSBR.dPz = zeros(numReconIso,1);
% Parameter derivatives: derivative w.r.t. R1 = 1/T1
solSBR.dR1x = zeros(numReconIso,1);
solSBR.dR1y = zeros(numReconIso,1);
solSBR.dR1z = zeros(numReconIso,1);
% Parameter derivatives: derivative w.r.t. R2 = 1/T2
solSBR.dR2x = zeros(numReconIso,1);
solSBR.dR2y = zeros(numReconIso,1);
solSBR.dR2z = zeros(numReconIso,1);
% Reference (acquired) signal: flatten array with order numRxCoils x numRxs
solSBR.Srefx = solFWD.Sx;
solSBR.Srefy = solFWD.Sy;
% Integrated signal: numRxCoils x numRxs
solSBR.Sx = 0*solFWD.Sx;
solSBR.Sy = 0*solFWD.Sy;
% Residual signal (Sref - S): numRxCoils x numRxs
solSBR.Rx = 0*solFWD.Sx;
solSBR.Ry = 0*solFWD.Sy;
% Gradients: numIsochromats  x numRxCoils
solSBR.GPr = zeros(numReconIso,numCoils);
solSBR.GPi = zeros(numReconIso,numCoils);
solSBR.GR1 = zeros(numReconIso,numCoils);
solSBR.GR2 = zeros(numReconIso,numCoils);
% Jacobians: numIsochromats x numCoils x numRxs
solSBR.JPDx = zeros(numReconIso,numCoils,numRxs);
solSBR.JPDy = zeros(numReconIso,numCoils,numRxs);
solSBR.JR1x = zeros(numReconIso,numCoils,numRxs);
solSBR.JR1y = zeros(numReconIso,numCoils,numRxs);
solSBR.JR2x = zeros(numReconIso,numCoils,numRxs);
solSBR.JR2y = zeros(numReconIso,numCoils,numRxs);


%% prepare function handle for the specific case
returnType  = 'jacobian'; % 'residual' /  'gradient' / 'jacobian'
pdActive    = 1;  % parameter active or not
r1Active    = 1;
r2Active    = 1;
fwdFun = @(x) spinTwin.sbrBloch.runFun( x, ...
    returnType, pdActive, r1Active, r2Active, ...
    solSBR, pulseSequence, sbrModel, simControl, dbgControl );

%% prepare initial solution
numIsochromats = sbrModel.numIsochromats;
numVars = (2*pdActive + r1Active + r2Active)*numIsochromats;
x = zeros(numVars,1);
numVars = 0;
if pdActive
    x(numVars+1:numVars+numIsochromats,1) = real(initialPD);
    numVars = numVars + numIsochromats;
    x(numVars+1:numVars+numIsochromats,1) = imag(initialPD);
    numVars = numVars + numIsochromats;
end
if r1Active
    x(numVars+1:numVars+numIsochromats,1) = 1/initialT1;
    numVars = numVars + numIsochromats;
end
if r2Active
    x(numVars+1:numVars+numIsochromats,1) = 1/initialT2;
    numVars = numVars + numIsochromats;
end

%% call the function
[fx, Gx] = fwdFun(x);

%% apply LS solution on PD (generates a Weighted image, same as FFT)
% crop Gx to remove the Jacobians of R1 and R2
Gx = Gx(:,1:2*numIsochromats);
% solve the LS problem to generate the Pr Pi solution
rSpace = Gx\fx;
% reorganize data: transform [Pr; Pi] into complex vector Pr+1j*Pi
rSpace = reshape(rSpace,[],2);
rSpace = rSpace(:,1) + 1j*rSpace(:,2);
% reshape into the recon sizes
rSpace = reshape(rSpace,reconNx,reconNy);

%% Compare to FFT 
% process the FFT image to remove th 2x in the X direction
% note that SBR recon was done on a 32x32 grid
figure();
subplot(1,2,1); imagesc(abs(iSpace(33:96,:))); title('FFT Image');
subplot(1,2,2); imagesc(abs(flipud(rSpace))); title('SBR Image');




