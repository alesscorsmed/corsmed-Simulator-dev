%
% SPINTWIN.TEST.EXAMPLE.DIFFUSION
%
%
%	Test diffusion effects on homogeneus cube using a PGSE sequence
%   90 -- +Gradient -- wait -- 180 -- +Gradient
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

clear all; close all; clc;

%% define parameters for the diffusion test
minDiffusion = 1e-9; maxDiffusion = 7e-9;
valuesTAU = [25, 100]*1e-3;
valuesAG  = [0.001, 0.030];
valuesTG  = [5, 20]*1e-3; 


%% initialize the sim Control with the default versioning for testing
[simControl] = spinTwin.setup.initializeSimControl();

%% prepare data
simControl.precision = 'single';
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

%% domain and model
b0 	 = 1.0;
mu 	 = 1.0;
fovX = 0.010;
fovY = 0.010;  
fovZ = 0.010;
dx   = 5e-4;
dy   = 5e-4;
dz   = 5e-4;
t1Value = 0.500;
t2Value = 0.250;
% generate spin Model
[spinModel] = spinTwin.test.model.homogeneousCubeModel( ...
    fovX, fovY, fovZ, dx, dy, dz, t1Value, t2Value, b0, mu );

simModel = spinModel.slice{1}.model;
%% get dimensions and initialize magnetizations
numCoils = simModel.numRxCoils;
numIso   = simModel.numIsochromats;

%% modify the diffusion values to be random for each voxel 
% within limits 1e-9 to 7e-9
xDiffusion = (1 + 6*rand(numIso,1))*1e-9;
yDiffusion = (1 + 6*rand(numIso,1))*1e-9;
zDiffusion = (1 + 6*rand(numIso,1))*1e-9;
% tissue type per voxel
tissueType = reshape(1:numIso,[],1);
% assign
simModel.numTissues         = numIso;
simModel.numProperties      = 6;
simModel.tissueValues       = zeros(numIso,6);
simModel.tissueValues(:,1)  = t1Value;
simModel.tissueValues(:,2)  = t2Value;
simModel.tissueType(:)      = reshape(tissueType,[],1);
% assign diffusion
simModel.tissueDiff      = zeros(numIso,3);
simModel.tissueDiff(:,1) = xDiffusion(:);
simModel.tissueDiff(:,2) = yDiffusion(:);
simModel.tissueDiff(:,3) = zDiffusion(:);

%% initial magnetizations
initialMagnitude = b0*mu;
initialPhase 	 = 0.0;
initialFA      	 = 0.0;
M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
M0p = initialPhase*ones(numIso,1);
M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
M0x = M0m.*cos(initialPhase);
M0y = M0m.*sin(initialPhase);

%% loop on PG parameters
testcase = 0;
for gradDir = ['x','y','z']
    for gradAmp = valuesAG
        for gradTime = valuesTG
            for TAU = valuesTAU
                
                %% update the values to generate the sequence
                grData.echoTime         = TAU;
                grData.duration      	= gradTime;
                grData.amplitude     	= gradAmp;
                grData.direction        = gradDir;
                
                %% get the sequence
                [ pulseSequence, grData ] = spinTwin.test.seq.pgse( ...
                    grData, rfData, mrSystem, dt);
                numRxs   = pulseSequence.numRxs;
                %% compute the reference ADC
                beta = grData.beta;
                switch lower(gradDir)
                    case 'x'
                        Dcoeff = simModel.tissueDiff(:,1);
                    case 'y'
                        Dcoeff = simModel.tissueDiff(:,2);
                    case 'z'
                        Dcoeff = simModel.tissueDiff(:,3);
                end
                ADCref = exp( -beta*abs(Dcoeff) );
                
                %% call the standard analytical kernel
                simControl.simulationEngine = 'Bloch';
                simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
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
                    solBloch, pulseSequence,simModel,...
                    motionModel, simControl, dbgControl );
                % magnitude and phase
                solBloch.Mm = abs(solBloch.Mx + 1j*solBloch.My);
                solBloch.Mp = angle(solBloch.Mx + 1j*solBloch.My);
                
                %% call the Phasor kernel (Analytical)
                simControl.simulationEngine = 'Diffusion';
                simControl.odeMethod        = 'explicit'; % analytical / explicit / implicit / adaptiveExp / adaptiveImp
                % initialize solution
                solDiffusion = data.simulation.initializeSolution(numIso,numRxs,numCoils);
                solDiffusion.indexes    = [];
                solDiffusion.Mx         = 1*M0x;
                solDiffusion.My         = 1*M0y;
                solDiffusion.Mz         = 1*M0z;
                solDiffusion.dMpDx      = 0*M0x;
                solDiffusion.dMpDy      = 0*M0x;
                solDiffusion.dMpDz      = 0*M0x;
                solDiffusion.Sx         = zeros(numCoils*numRxs,1);
                solDiffusion.Sy         = zeros(numCoils*numRxs,1);
                solDiffusion.Sz         = zeros(numCoils*numRxs,1);
                % kernel call
                [solDiffusion] = spinTwin.fwdBloch.runDiffusion(...
                    solDiffusion, pulseSequence,simModel,...
                    motionModel, simControl, dbgControl );
                % magnitude and phase
                solDiffusion.Mm = abs(solDiffusion.Mx + 1j*solDiffusion.My);
                solDiffusion.Mp = angle(solDiffusion.Mx + 1j*solDiffusion.My);
                
                %% compute numerical ADC
                ADCsim = solDiffusion.Mm./solBloch.Mm;
                
                %% plot data
                testcase = testcase + numIso;
                Bcoeff = beta*ones(size(ADCsim));
                figure(1); hold on;
                plot(Bcoeff(:), ADCref(:), 'go', 'LineWidth', 2, 'MarkerSize',5 );
                plot(Bcoeff(:), ADCsim(:), 'b+', 'LineWidth', 1,'MarkerSize',5 );
                LogErr = log10(abs(ADCsim-ADCref));
                LogErr(LogErr < -6) = -6;
                plot(Bcoeff(:), LogErr(:), 'r.', 'MarkerSize',10 );
                xlabel('B coefficient'); ylabel('ADC'); grid on;
                ylim([-6.0001,1.0001]);
                legend('Theoretical', 'Simulated', 'Abs. Error (Logscale)')
                title(sprintf('Diffusion test: %d cases',testcase));
                
            end
        end
    end
end

