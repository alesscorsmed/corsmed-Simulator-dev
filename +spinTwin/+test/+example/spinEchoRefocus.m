%
% SPINTWIN.TEST.EXAMPLE.SPINECHOREFOCUS
%
%	Script to verify the effect of a 180 refocusing pulse 
%   when initial phase (and its derivative) is not 0.
%   Generates a SE encoding (90 + GR + 180 + GR), and
%   checks the echo signal.
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
simControl.precision = 'single';

% gradient and echo time
phaseEncodings              = 0; % array with encodings for repetitions, e.g. linspace(-1,1,11);
acqData.preEncoding         = 1; % if 1 applies FE RW before 180, if 0, applies FE RW after 180
acqData.echoTime            = 15e-3;
acqData.repTime             = 4.000; % TR
acqData.rxBW                = 50e3;  % in Hz
acqData.fovFE               = 0.3; % set these to 0.05 to see some spurious echoes
acqData.fovPE               = 0.3;
acqData.fovSE               = 0.3;
acqData.numFE               = 128;
acqData.numPE               = 128;
acqData.numSE               = 128;
acqData.samplingFactorFE    = 2;
acqData.samplingFactorPE    = 1;
acqData.samplingFactorSE    = 1;
% RF (180 refocusing pulse is the same, but with 180 FA)
rfData.type             = 'sinc';
rfData.flipAngle      	= 90; % in degrees
rfData.phase         	= 0.0; % in rads
rfData.duration         = 3e-3;
rfData.cycles           = 2;
rfData.sliceThickness   = 6e-3;
rfData.doSliceSelect 	= 0; % set this to 1 for slice selection
rfData.preRWScale    	= 0; % set this to 1 for slice selection
rfData.postRWScale    	= 0; % set this to 1 for slice selection
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
[ pulseSequence ] = spinTwin.test.seq.se( ...
    acqData, rfData, mrSystem, dt, phaseEncodings );

%% domain and model
b0 	 = 1.0;
mu 	 = 1.0;
fovX = 0.050;
fovY = 0.050;  
fovZ = 0.010;
dx   = 1e-3;
dy   = 1e-3;
dz   = 1e-3;
t1Limits = [0.500, 0.500]; % put a range if you want to have various T1/T2
t2Limits = [0.100, 0.100];
% generate spin Model
[spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
    fovX, fovY, fovZ, dx, dy, dz, mean(t1Limits), mean(t2Limits), b0, mu );

% get and modify the simulation model
simModel = spinModel.slice{1}.model;
if ( numel(unique(t1Limits)) > 1 ) && ( numel(unique(t2Limits)) > 1 )
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
end

%% get dimensions and initialize magnetizations
numCoils = simModel.numRxCoils;
numIso   = simModel.numIsochromats;
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
    solBloch, pulseSequence, simModel,...
    motionModel, simControl, dbgControl );
% magnitude and phase
solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);


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
    solPhasor, pulseSequence, simModel,...
    motionModel, simControl, dbgControl );

% magnitude and phase
solPhasor.Mm = abs(solPhasor.Mx + 1j*solPhasor.My);
solPhasor.Mp = angle(solPhasor.Mx + 1j*solPhasor.My);


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
    solPhasorODE, pulseSequence, simModel,...
    motionModel, simControl, dbgControl );

% magnitude and phase
solPhasorODE.Mm = abs(solPhasorODE.Mx + 1j*solPhasorODE.My);
solPhasorODE.Mp = angle(solPhasorODE.Mx + 1j*solPhasorODE.My);



%% plot sequence
figure();
plot(pulseSequence.time,pulseSequence.rfmSignal*1e3);
xlabel('time (s)');
hold on
plot(pulseSequence.time,pulseSequence.gxSignal);
plot(pulseSequence.time,pulseSequence.gySignal);
plot(pulseSequence.time,pulseSequence.gzSignal);
plot(pulseSequence.time(pulseSequence.rxSignal>0), ...
    pulseSequence.gxSignal(pulseSequence.rxSignal>0), 'o');
plot(pulseSequence.time(pulseSequence.rxLimits(:,1)), ...
    pulseSequence.gxSignal(pulseSequence.rxLimits(:,1)), '>');
plot(pulseSequence.time(pulseSequence.rxLimits(:,2)), ...
    pulseSequence.gxSignal(pulseSequence.rxLimits(:,2)), '<');
plot(pulseSequence.time(pulseSequence.partLimits(:,1)),zeros(pulseSequence.numParts,1), '^');
plot(pulseSequence.time(pulseSequence.partLimits(:,2)),zeros(pulseSequence.numParts,1), 'v');
if nnz(pulseSequence.swcSignal) > 0
    plot(pulseSequence.time(pulseSequence.swcSignal>0), ...
        zeros(nnz(pulseSequence.swcSignal),1), 's', 'LineWidth', 2, 'MarkerSize', 10);
    legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end', 'SWC');
else
    legend('RF amp', 'Gx', 'Gy', 'Gz', 'RXs', 'RX start', 'Rx end', 'Part Start', 'Part end');
end


%% plot signal
figure();
scale = max(abs(solPhasorODE.Sx + 1j*solPhasorODE.Sy));
rxTimePoints = pulseSequence.time(pulseSequence.rxSignal>0);

figAx(1) = subplot(3,1,1);
plot(rxTimePoints,solBloch.Sx(:));
hold on
plot(rxTimePoints,solPhasor.Sx(:));
plot(rxTimePoints,solPhasorODE.Sx(:));
xlabel('time (s)'); ylabel('Sx');
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
axis([min(rxTimePoints), max(rxTimePoints), -scale, +scale]);

figAx(2) = subplot(3,1,2);
plot(rxTimePoints,solBloch.Sy(:));
hold on
plot(rxTimePoints,solPhasor.Sy(:));
plot(rxTimePoints,solPhasorODE.Sy(:));
xlabel('time (s)'); ylabel('Sy');
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
axis([min(rxTimePoints), max(rxTimePoints), -scale, +scale]);

figAx(3) = subplot(3,1,3);
plot(rxTimePoints,abs(solBloch.Sx + 1j*solBloch.Sy));
hold on
plot(rxTimePoints,abs(solPhasor.Sx + 1j*solPhasor.Sy));
plot(rxTimePoints,abs(solPhasorODE.Sx + 1j*solPhasorODE.Sy));
xlabel('time (s)'); ylabel('|S|');
legend('Analytic','Phasor (Analytic)', 'Phasor (ODE)');
axis([min(rxTimePoints), max(rxTimePoints), -scale, +scale]);

linkaxes(figAx, 'x');

