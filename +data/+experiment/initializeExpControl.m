function [expControl] = initializeExpControl(approach)
%
% DATA.EXPERIMENT.INITIALIZE
%
%	Function that initializes the expControl data structure,
%   with information for the experiment control.
%   Returns a simulationControl with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   sessionData   struct with sessionData
%
% OUTPUT
%   expControl   expControl structure
%
%========================  CORSMED AB © 2020 ==============================
%

if nargin<1
    approach = '';
end

%% default data
expControl.name         = 'experiment';
expControl.userID       = '';
expControl.instanceID   = '';
expControl.experimentID = '';
expControl.courseID     = '';
expControl.pulseqID     = '';
expControl.versionNum   = '';
expControl.approach     = approach;

%% time stamp
timeStamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
expControl.timeStamp    = timeStamp; % display format for messages
% change time stamp into file format yyyymmddHHMMSS
timeStamp = strrep(timeStamp,'-','');
timeStamp = strrep(timeStamp,':','');
timeStamp = strrep(timeStamp,' ','');
expControl.fileTimeStamp = timeStamp; % file format for create files and ids

%% connections and folder systems
expControl.progress     = 0;
if strcmp(approach,'standalone')
    expControl.connLocalDB  = [];
end
expControl.folderSystem = [];

%% debugging and development mode
expControl.debug.devMode      = 1;
expControl.debug.debugMode    = 1;
expControl.debug.debugFile    = [];
expControl.debug.waitBarBE    = 0;

%% model related
expControl.model.useSliceThickness  = 1;
expControl.model.pdInhomogeneity    = 1;
expControl.model.isotropicGrid      = 1;
expControl.model.coilType           = 'optimal';
expControl.model.gridSizeOption     = '1mm';
expControl.model.gridStep           = [1.0, 1.0, 1.0]*1e-3;
expControl.model.coilMode           = 'basic';  % It refers to the option coil_basic_or_advanced
expControl.model.xExtendPct         = 0.1; % extend percentage in x in model interpolation (10% default)
expControl.model.rxCoilName         = 'none'; % no coil selected yet
expControl.model.zIsocenter         = 0.0; % isocenter for the coil positioning
expControl.model.cardiacPhase       = 1; % index of the cardiac phase to use

%% perfusion
% parameters
expControl.model.perfusion.GdRelaxivity	= 5.6 / 1000;   % Gd relaxivity         [l/(mmol*ms)]
expControl.model.perfusion.MBFrest      = 1.0 / 60;     % Rest MBF              [ml/g/s]
expControl.model.perfusion.MBFstress    = 3.5 / 60;     % Stress MBF            [ml/g/s]
expControl.model.perfusion.Fermi_alpha  = 0.25;         % Fermi model parameter alpha
expControl.model.perfusion.Fermi_beta   = 0.25;         % Fermi model parameter beta
expControl.model.perfusion.Tshift       = 3;            % temporal LV-myo shift [s]
expControl.model.perfusion.contrastDose = 0.075;        % reference dose [mmol/kg]
expControl.model.perfusion.RestStress   = 2;            % 1=rest; 2=stress
% perfursion tissues ID
expControl.model.perfusion.idMyoLV      = 1500;
expControl.model.perfusion.idMyoLA      = 1497;
expControl.model.perfusion.idMyoRV      = 1499;
expControl.model.perfusion.idMyoRA      = 1498;
expControl.model.perfusion.idBloodLV    = 1504;
expControl.model.perfusion.idBloodLA    = 1501;
expControl.model.perfusion.idBloodRV    = 1503;
expControl.model.perfusion.idBloodRA    = 1502;
%  ids of other tissues to keep intact (non zero)
expControl.model.perfusion.idSteady     = [40];
% timing     
expControl.model.perfusion.contrasts    = 32;  % Number of Dynamics
expControl.model.perfusion.Tcc          = 1.0; % Cardiac cycle duration [s]
% flag to apply
expControl.model.perfusion.apply        = 0;


%% simulation related
expControl.simulation.numberOfSim       = 0; % number of simulations
expControl.simulation.analyticalSim     = 0; % full analytical
expControl.simulation.simulationEngine  = 'analytical'; % analytical / phasor / diffusion / numerical
expControl.simulation.precision         = 'single'; % single / double
expControl.simulation.odeMethod         = 'adaptiveExp'; % explicit / adaptiveExp / implicit / adaptiveImp
expControl.simulation.simulationKernel  = 'latest'; % latest / stable
expControl.simulation.activateCS        = 0;  % activate chemical shift
expControl.simulation.activateSusc      = 0;  % activate susceptibility
expControl.simulation.applyNoise        = 0;

%% parallelization and GPU mode
expControl.simulation.parEval       = 0;
expControl.simulation.numGPUs    	= 1;
expControl.simulation.threads     	= 256;
expControl.simulation.blocks     	= 256;
expControl.simulation.kernelVer  	= 'v26_20210128_ufm';
expControl.simulation.architecture  = 'sm70';
expControl.simulation.kernelFolder  = '/efs-mount-point/S20/PROJECTS/edutool/kernels/';
expControl.simulation.kernelPtx  	= sprintf('%s%s_%s.ptx', ...
    expControl.simulation.kernelFolder, expControl.simulation.kernelVer,...
    expControl.simulation.architecture);

%% sequence related
expControl.sequence.dtGR            = 10e-6; % default time discretization for GR
expControl.sequence.dtRF            = 5e-6; % default time discretization for RF
expControl.sequence.tGuardRF        = 10e-6; % time guard to avoid RF/GR overlap
expControl.sequence.dwellTime       = 'dynamic';
expControl.sequence.timeCompression = 0; % former fast algorithm
expControl.sequence.minContextExc   = 0; % minimizes context change in kernel exec
expControl.sequence.deactivateGx    = 0;
expControl.sequence.deactivateGy    = 0;
expControl.sequence.deactivateGz    = 0;
expControl.sequence.deactivateSS    = 0; % deactivate slice selection

%% motion related
expControl.motionSpecs.pattern      = 0;
expControl.motionSpecs.rotAngle     = 0;
expControl.motionSpecs.rotFreq      = 0;
expControl.motionSpecs.transMag     = 0;
expControl.motionSpecs.transFreq    = 0;
expControl.motionSpecs.transAxis    = 0;

%% standard and other
expControl.advancedNotifications   = 0;
expControl.b0map                   = 0;
%% system
expControl.mrSystem = data.experiment.initializeMRSystem();
