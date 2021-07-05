%% initial setup
% json filename
jsonPath        = ['+sbr',filesep,'testInputs',filesep,'testInput_fromNewSimulator_20210122a_v26_sm70.json'];
% jsonPath        = ['inputs/testInput_20210119_singleGradient_singlePrecision.json'];
% jsonPath = ['/efs-mount-point-MATLAB/cxanthis/compiled_projects/sbrMain/for_redistribution_files_only/',...
%     'inputs/testInput_20210122_singlePrecision_doubleGradient_v26_sm70.json'];

% ground-truth kspace path
kSpacePath      = ['inputs',filesep,'kspaceGT_fromNewSimulator.mat'];

% anatomical model path
anatModelPath   = '';

% pulse sequence path
pulseSeqPath    = ['/efs-mount-point-MATLAB/cxanthis/corsmed-Simulator/+sbr',filesep,'testInputs',filesep,...
    'bSSFP_128x128_7_6_sineFA_k_10_singleLobe_withoutaby2_200kHz_CXedit.mat'];

%% Initial estimation
x0         = [1.26 0.21 0.95 1.26 0.21 0.95 1.26 0.21 0.95 1.26 0.21 0.95];     % The starting point.
options.lb = [1.2 0.2 0.9 1.2 0.2 0.9 1.2 0.2 0.9 1.2 0.2 0.9];                 % Lower bound on the variables.
options.ub = [1.8 0.3 1.1 1.8 0.3 1.1 1.8 0.3 1.1 1.8 0.3 1.1];                 % Upper bound on the variables.

%% Set the IPOPT options.  
options.ipopt.hessian_approximation         = 'limited-memory';
options.ipopt.limited_memory_update_type    = 'bfgs';
options.ipopt.derivative_test               = 'none';
options.ipopt.linear_solver                 = 'mumps'; %'ma57'; 'pardiso';
options.ipopt.mu_strategy                   = 'adaptive';
options.ipopt.print_level                   = 5;
options.ipopt.tol                           = 1e-3;
options.ipopt.max_iter                      = 20;

% additional costumization
options.ipopt.mumps_pivtol                  = 0.0001;

%% Set the specs for the gradient function
options.auxdata.gradient.stepT1 = 0.04;
options.auxdata.gradient.stepT2 = 0.002;
options.auxdata.gradient.stepPD = 0.02;
options.auxdata.gradient.side   = 'double'; % single / double % refers to single-side or double-side gradient


%% define the characteristics of the pulse sequence
% If pulseSequence.type is empty, the pulse sequence will be loaded from
% the pulseSequence.path. If not, the specific pulse sequence type will
% be created based on the pulseSequence.specs.
options.auxdata.pulseSequence.type          = '';
options.auxdata.pulseSequence.specs         = '';
options.auxdata.pulseSequence.path          = pulseSeqPath;


%% Run the ground-truth simulation or load the ground-truth kspace
% if runSimulationGT = 1, then run the GT simulation with the 
% options.auxdata.simulationGT.anatomicalPath and the 
% options.auxdata.pulseSequence. Otherwise, load the ground-truth kspace
options.auxdata.runSimulationGT = 1;

% Run the ground-truth experiment. Two options: load the anatomical model
% from path or create a new one
options.auxdata.simulationGT.loadAnatModel                      = 0;
options.auxdata.simulationGT.anatomicalPath                     = '';

options.auxdata.simulationGT.anatomicalModel.specs.xMin         = -0.2;   % in m
options.auxdata.simulationGT.anatomicalModel.specs.xMax         = 0.2;    % in m
options.auxdata.simulationGT.anatomicalModel.specs.yMin         = -0.2;   % in m
options.auxdata.simulationGT.anatomicalModel.specs.yMax         = 0.2;    % in m
options.auxdata.simulationGT.anatomicalModel.specs.zMin         = 0;      % in m
options.auxdata.simulationGT.anatomicalModel.specs.zMax         = 0;      % in m
options.auxdata.simulationGT.anatomicalModel.specs.gridStepX    = 0.001;  % in m
options.auxdata.simulationGT.anatomicalModel.specs.gridStepY    = 0.001;  % in m
options.auxdata.simulationGT.anatomicalModel.specs.gridStepZ    = 0.001;  % in m
options.auxdata.simulationGT.anatomicalModel.specs.tissuesX     = 2;
options.auxdata.simulationGT.anatomicalModel.specs.tissuesY     = 2;
options.auxdata.simulationGT.anatomicalModel.specs.tissueValues = ...
    [1.5,0.25,1,0,0,0;...
    1.5,0.25,1,0,0,0;...
    1.5,0.25,1,0,0,0;...
    1.5,0.25,1,0,0,0];

% define the path to the ground-truth kspace
% The kspace should be saved in a certain matrix-name: 'kSpace'
options.auxdata.simulationGT.kSpace.path                        = kSpacePath;


%% define the characteristics of the anatomical model
% If anatomicalModel.path is empty, the anatomical model will be created 
% based on the anatomicalModel.specs. If not, the anatomical model will be
% loaded.
options.auxdata.anatomicalModel.path            = anatModelPath;
options.auxdata.anatomicalModel.specs.xMin      = -0.2;   % in m
options.auxdata.anatomicalModel.specs.xMax      = 0.2;    % in m
options.auxdata.anatomicalModel.specs.yMin      = -0.2;   % in m
options.auxdata.anatomicalModel.specs.yMax      = 0.2;    % in m
options.auxdata.anatomicalModel.specs.zMin      = 0;      % in m
options.auxdata.anatomicalModel.specs.zMax      = 0;      % in m
options.auxdata.anatomicalModel.specs.gridStepX = 0.001;  % in m
options.auxdata.anatomicalModel.specs.gridStepY = 0.001;  % in m
options.auxdata.anatomicalModel.specs.gridStepZ = 0.001;  % in m
options.auxdata.anatomicalModel.specs.tissuesX  = 2;
options.auxdata.anatomicalModel.specs.tissuesY  = 2;

%% define the characteristics of the experiment
options.auxdata.sessionData.errorFolder         = 'errors/';            
options.auxdata.sessionData.kernelFolder        = 'kernels/';
options.auxdata.sessionData.kernelPtx           = '/efs-mount-point-MATLAB/cxanthis/compiled_projects/sbrMain/for_redistribution_files_only/kernels/v26_sm70.ptx';
options.auxdata.sessionData.simulationEngine    = 'analytical'; % analytical / phasor / diffusion / numerical
options.auxdata.sessionData.threads             = 256;
options.auxdata.sessionData.blocks              = 256;
options.auxdata.sessionData.debugMode           = 0;
options.auxdata.sessionData.debugFile           = '';
options.auxdata.sessionData.numGPU              = 1;
options.auxdata.sessionData.precision           = 'single'; % single / double

%% Save to json file
experimentData.options                      = options;
experimentData.x0                           = x0;

jsonStr                                     = jsonencode(experimentData);

fid = fopen(jsonPath, 'w');
if fid == -1
    error('Cannot create JSON file'); 
end

fwrite(fid, jsonStr, 'char');
fclose(fid);