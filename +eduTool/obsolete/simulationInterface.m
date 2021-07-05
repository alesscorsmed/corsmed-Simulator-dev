function [Sx,Sy,Sz] = simulationInterface(...
    anatModel_struct,pulseSeq_struct,expControl)
%
%  Interface to nSim.py, that runs the numerical simulation.
%  calls directly the python function passing arrays
%

% %% control for the experiment
% [expControl] = data.expControl.initialize();
% expControl.kernelPtx = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_23.ptx';

%% generate spin model and populate
fprintf(1, 'Preparing anatomical model ... ');
t_start = tic;
% 1 slice
[spinModel] = data.spinModel.initialize(1);
niso = size(anatModel_struct.model_spatial,1);
% assign data
spinModel.name                          = anatModel_struct.modelName;
spinModel.totalIsochromats              = niso;

spinModel.slice{1}.model.resolution     = [anatModel_struct.resolution(1);...
    anatModel_struct.resolution(2); anatModel_struct.resolution(3)];
spinModel.slice{1}.model.numIsochromats = niso;
spinModel.slice{1}.model.nonZeroIndex   = 1:niso;

spinModel.slice{1}.model.r3D            = anatModel_struct.model_spatial;
spinModel.slice{1}.model.x              = anatModel_struct.model_spatial(:,1);
spinModel.slice{1}.model.y              = anatModel_struct.model_spatial(:,2);
spinModel.slice{1}.model.z              = anatModel_struct.model_spatial(:,3);

spinModel.slice{1}.model.b0             = 1.0;
spinModel.slice{1}.model.mu             = 1.0;
if isfield(anatModel_struct, 'inhom') ...
        && (length(anatModel_struct.inhom(:)) == niso)
    spinModel.slice{1}.model.bi = single(reshape(anatModel_struct.inhom,niso,1));
else
    spinModel.slice{1}.model.bi = single(zeros(niso,1));
end
if isfield(anatModel_struct, 'PDinhom') ...
        && (length(anatModel_struct.PDinhom(:)) == niso)
    spinModel.slice{1}.model.pd = single(reshape(anatModel_struct.PDinhom,niso,1));
else
    spinModel.slice{1}.model.pd = single(ones(niso,1));
end

spinModel.slice{1}.model.xDiffusion     = zeros(niso,1);
spinModel.slice{1}.model.yDiffusion     = zeros(niso,1);
spinModel.slice{1}.model.zDiffusion     = zeros(niso,1);

spinModel.slice{1}.model.xDiffusion(anatModel_struct.model_new == 12) = 10e-9;

spinModel.slice{1}.model.numTissues     = size(anatModel_struct.model_tissues,1); 
spinModel.slice{1}.model.numProperties  = size(anatModel_struct.model_tissues,2); 
spinModel.slice{1}.model.tissueValues   = anatModel_struct.model_tissues.';
spinModel.slice{1}.model.tissueType     = anatModel_struct.model_new;

% coils
if isfield(anatModel_struct, 'coilmapsx')
    spinModel.slice{1}.model.coilMapsX = single(anatModel_struct.coilmapsx);
else
    spinModel.slice{1}.model.coilMapsX = single(ones(niso,1));
end
if isfield(anatModel_struct, 'coilmapsy')
    spinModel.slice{1}.model.coilMapsY = single(anatModel_struct.coilmapsy);
else
    spinModel.slice{1}.model.coilMapsY = single(zeros(niso,1));
end
spinModel.slice{1}.model.numCoils  = size(spinModel.slice{1}.model.coilMapsX,2);

% and done
fprintf(1, ' done -- elapsed time %f\n', toc(t_start));

%% Preparing the sequence data
fprintf(1, 'Processing sequence ... ');
t_start = tic;

if isempty(pulseSeq_struct.newPulseSequence)
    % initialize
    pulseSequence = data.pulseSequence.initialize();
    
    % basic info
    pulseSequence.name       = 'N/A';
    pulseSequence.type       = 'N/A';
    pulseSequence.endEvent   = 'none'; % indicates what happens at the end
    pulseSequence.gamma      = 42.577478518e6;
    pulseSequence.numSteps   = size(pulseSeq_struct.pulse_sequence,2);
    
    pulseSequence.time      = zeros(pulseSequence.numSteps,1); % times
    pulseSequence.timeDiff  = zeros(pulseSequence.numSteps,1); % time deltas
    pulseSequence.rxSignal  = zeros(pulseSequence.numSteps,1); % receiver readout
    pulseSequence.swcSignal = zeros(pulseSequence.numSteps,1); % software crusher
    pulseSequence.gxSignal  = zeros(pulseSequence.numSteps,1); % x gradient
    pulseSequence.gySignal  = zeros(pulseSequence.numSteps,1); % y gradient
    pulseSequence.gzSignal  = zeros(pulseSequence.numSteps,1); % z gradient
    pulseSequence.rfmSignal = zeros(pulseSequence.numSteps,1); % RF magnitude
    pulseSequence.rfpSignal = zeros(pulseSequence.numSteps,1); % RF phase
    pulseSequence.rffSignal = zeros(pulseSequence.numSteps,1); % RF frequency
    
    pulseSequence.timeDiff(:)   = (pulseSeq_struct.pulse_sequence(8,:)-...
        pulseSeq_struct.pulse_sequence(7,:)+1)*pulseSeq_struct.dt;
    pulseSequence.time(:)       = cumsum(pulseSequence.timeDiff);
    
    pulseSequence.rxSignal(:)   = pulseSeq_struct.isInKspace(:);
    pulseSequence.swcSignal(:)  = pulseSeq_struct.soft_crushers(:);
    
    pulseSequence.rfmSignal(:)  = pulseSeq_struct.pulse_sequence(1,:);
    pulseSequence.rfpSignal(:)	= pulseSeq_struct.pulse_sequence(2,:);
    pulseSequence.rffSignal(:)  = pulseSeq_struct.pulse_sequence(3,:);
    pulseSequence.gxSignal(:)   = pulseSeq_struct.pulse_sequence(4,:);
    pulseSequence.gySignal(:)	= pulseSeq_struct.pulse_sequence(5,:);
    pulseSequence.gzSignal(:)   = pulseSeq_struct.pulse_sequence(6,:);
    
    pulseSequence.totalTime     = pulseSequence.time(end);
    pulseSequence.numRxs        = nnz(pulseSequence.rxSignal(:));
    
    pulseSequence.numParts     = 1; % number of parts
    pulseSequence.partType{1}  = 'DF'; % type of part: RF / RO / GR / DF
    pulseSequence.partLimits   = [1, pulseSequence.numSteps]; % index start/end of part
    
    expControl.engine          = 'analytical';
    
else
    pulseSequence = pulseSeq_struct.newPulseSequence;
end

% and done
fprintf(1, ' done -- elapsed time %f\n', toc(t_start));

%% initialize the solution
fprintf(1, 'Initializing solution ... ');
t_start = tic;
solution = data.solution.initialize(spinModel.slice{1}.model.numIsochromats,...
    pulseSequence.numRxs,spinModel.slice{1}.model.numCoils);
% initialize initial conditions
solution.Mz(:) = spinModel.slice{1}.model.b0*spinModel.slice{1}.model.mu*...
    spinModel.slice{1}.model.pd(:);
% and done
fprintf(1, ' done -- elapsed time %f\n', toc(t_start));

%% call the simulator
switch expControl.engine
    case 'phasor'
        [solution] = simulator.bloch.kernel.runPhasor(solution,...
            pulseSequence,spinModel.slice{1}.model,expControl);
    case 'diffusion'
        [solution] = simulator.diffusion.kernel.runSDW(solution,...
            pulseSequence,spinModel.slice{1}.model,expControl);
    otherwise % 'analytical'
        [solution] = simulator.bloch.kernel.runAnalytical(solution,...
            pulseSequence,spinModel.slice{1}.model,expControl);
end

%% process solution
Sx = reshape(solution.Sx,spinModel.slice{1}.model.numCoils,pulseSequence.numRxs).';
Sy = reshape(solution.Sy,spinModel.slice{1}.model.numCoils,pulseSequence.numRxs).';
Sz = reshape(solution.Sz,spinModel.slice{1}.model.numCoils,pulseSequence.numRxs).';

