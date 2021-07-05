%
% SIMULATOR.EXAMPLE.DIFFUSIONPGSEEXAMPLE
%
%	Test diffusion effects on homogeneus cube using a PGSE sequence
%   90 -- +Gradient -- wait -- 180 -- +Gradient
%
% INPUT
%   see below for the parameters to play with
%
%========================  CORSMED AB © 2020 ==============================
%

%% define parameters for the diffusion test
minDiffusion = 1e-9; maxDiffusion = 7e-9;
valuesTAU = [25, 100]*1e-3;
valuesAG  = [0.001, 0.030];
valuesTG  = [5, 20]*1e-3; 

%% define RF conditions
gamma           = 42.577478518e6;
b0              = 1.0;
mu              = 1;
f0              = gamma*b0;
dt              = 1e-6;
maxGStrenght    = 40e-3;
gradSlewrate    = 150.0;
sliceThickness  = 6e-3;
rfDuration      = 3e-3;
rfFlipAngle     = 90;
rfPhase         = 0.0;
rfCycles        = 2;
preRWScale      = 0;
postRWScale     = 0;
doSliceSelect   = 0;

%% initial magnetization state
initialMagnitude    = b0*mu;
initialPhase        = 0.0;
initialFA           = 0.0;

%% domain
fovX = 0.010;
fovY = 0.010;  
fovZ = 0.010;
dx   = 5e-4;
dy   = 5e-4;
dz   = 5e-4;
t1Value = 0.500;
t2Value = 0.250;

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

%% initial magnetizations
% get sizes
numCoils = spinModel.slice{1}.model.numCoils;
numIso   = spinModel.slice{1}.model.numIsochromats;
numRxs   = 1;

M0m = initialMagnitude*sin(initialFA*pi/180)*ones(numIso,1);
M0p = initialPhase*ones(numIso,1);
M0z = initialMagnitude*cos(initialFA*pi/180)*ones(numIso,1);
M0x = M0m.*cos(initialPhase);
M0y = M0m.*sin(initialPhase);

%% modify the diffusion values to be random for each voxel 
% within limits 1e-9 to 7e-9
spinModel.slice{1}.model.xDiffusion(:) = (1 + 6*rand(numIso,1))*1e-9;
spinModel.slice{1}.model.yDiffusion(:) = (1 + 6*rand(numIso,1))*1e-9;
spinModel.slice{1}.model.zDiffusion(:) = (1 + 6*rand(numIso,1))*1e-9;


%% loop on PG parameters
testcase = 0;
for gradDir = ['x','y','z']
    for gradAmp = valuesAG
        for gradTime = valuesTG
            for TAU = valuesTAU
                      
                %% assign values
                pulseSequence = data.pulseSequence.initialize();
                pulseSequence.encPG.Dir = gradDir;
                pulseSequence.encPG.TG  = gradTime;
                pulseSequence.encPG.AG  = gradAmp;
                pulseSequence.encPG.TAU = TAU;
                % main RF info
                pulseSequence.mainRF.flipAngle          = 90; % in degrees
                pulseSequence.mainRF.phase              = -pi/2; % in rads
                pulseSequence.mainRF.duration           = rfDuration;
                pulseSequence.mainRF.cycles             = rfCycles;
                pulseSequence.mainRF.sliceThickness     = sliceThickness;
                pulseSequence.mainRF.doSliceSelect      = doSliceSelect;
                pulseSequence.mainRF.preRWScale         = preRWScale;
                pulseSequence.mainRF.postRWScale        = postRWScale;
                % refocusing RF info
                pulseSequence.refRF.flipAngle           = 180; % in degrees
                pulseSequence.refRF.phase               = 0; % in rads
                pulseSequence.refRF.duration            = rfDuration;
                pulseSequence.refRF.cycles              = rfCycles;
                pulseSequence.refRF.sliceThickness      = sliceThickness;
                pulseSequence.refRF.doSliceSelect       = doSliceSelect;
                pulseSequence.refRF.preRWScale          = preRWScale;
                pulseSequence.refRF.postRWScale         = postRWScale;
                
                %% generate the sequence
                [pulseSequence] = sequence.pgse.encodingPGSE(...
                    pulseSequence, mrSystem, expControl);
                
                % cheat: modify sequence to have 1 RO only at the end
                pulseSequence.numRxs                  = 1;
                pulseSequence.subSeq{1}.numRxs        = 1;
                pulseSequence.subSeq{1}.rxSignal(:)   = 0;
                pulseSequence.subSeq{1}.rxSignal(end) = 1;
                
                %% compute the reference ADC
                beta = pulseSequence.encPG.Beta;
                switch lower(gradDir)
                    case 'x'
                        Dcoeff = spinModel.slice{1}.model.xDiffusion;
                    case 'y'
                        Dcoeff = spinModel.slice{1}.model.yDiffusion;
                    case 'z'
                        Dcoeff = spinModel.slice{1}.model.zDiffusion;
                end
                ADCref = exp( -beta*abs(Dcoeff) );
                
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
                
                %% call the diffusion kernel
                solDiff = data.solution.initialize(numIso,numRxs,numCoils);
                solDiff.indexes    = [];
                solDiff.Mm         = 1*M0m;
                solDiff.Mp         = 1*M0p;
                solDiff.Mx         = 1*M0x;
                solDiff.My         = 1*M0y;
                solDiff.Mz         = 1*M0z;
                solDiff.dMpDx      = 0*M0p;
                solDiff.dMpDy      = 0*M0p;
                solDiff.dMpDz      = 0*M0p;
                solDiff.Sx         = zeros(numCoils*numRxs,1);
                solDiff.Sy         = zeros(numCoils*numRxs,1);
                solDiff.Sz         = zeros(numCoils*numRxs,1);
                % kernel call
                [solDiff] = simulator.diffusion.kernel.runSDW(solDiff,...
                    pulseSequence,spinModel.slice{1}.model,expControl);
                
                %% compute numerical ADC
                ADCsim = solDiff.Mm./solAnalytical.Mm;
                
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

