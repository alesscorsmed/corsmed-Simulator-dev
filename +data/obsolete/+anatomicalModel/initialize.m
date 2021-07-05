function [anatomicalModel] = initialize(anatomicalID,basePath,application)
%
% DATA.ANATOMICALMODEL.INITIALIZE
%
%	Function that initializes an anatomical model.
%   Loads the data from file.
%   If not file is passed, returns a model with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   filename   path of the data to load in the model
%
% OUTPUT
%   anatomicalModel   anatomicalModel structure
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'data.anatomicalModel.initialize';
if (nargin < 1)
    anatomicalID = 6;
end
if (nargin < 2) || isempty(basePath)
    ME = MException('Model:wrongPath',...
        '%s : Anatomical model path not available',functionName);
    throw(ME);
end
if (nargin < 3) || isempty(application)
    application = 'sbr';
end


%% report start
fid = 1;
tTotal = tic();
fprintf(fid, '\n%s : start', functionName);
fprintf(fid, '\n  Loading data from %s', basePath);

%% load the data
anatomicalModel = [];
% % % data inside
% % anatomicalModel.name                = 'MIDA';
% % anatomicalModel.b0                  = 1.0;
% % anatomicalModel.mu                  = 1.0;
% % anatomicalModel.numVoxels           = 80640000; 
% % anatomicalModel.backgroundTissue    = 50;
% % anatomicalModel.fatTissuesIDs       = [43,62];
% % anatomicalModel.resolution          = [0.5, 0.5, 0.5]*1e-3;
% % anatomicalModel.dimensions          = [480, 350, 480];
% % anatomicalModel.domain              = [0.1750 0.2423 0.2400];
% % anatomicalModel.spatial             = []; % numVoxels x 1
% % anatomicalModel.numTissues          = 116;
% % anatomicalModel.numProperties       = 6;
% % anatomicalModel.tissueType          = []; % numVoxels x 1
% % anatomicalModel.tissueValues        = []; % numTissues x numProperties
% % anatomicalModel.tissueDiff          = []; % numTissues x 3
% % anatomicalModel.pdFluctuationFactor = 0.05;
% % anatomicalModel.pdInhomogeneity     = []; % numVoxels x 1
% % anatomicalModel.b0Inhomogeneity     = []; % numVoxels x 1
% % anatomicalModel.isGridded           = 1;

switch anatomicalID
    
    case 4 % XCAT model
        if strcmpi(application, 'edutool')
            fileName = 'XCAT_20210110.mat';
            load( sprintf('%s%s',basePath,fileName), 'anatomicalModel' );
        else
            fileName = 'anatomicalModel_1mm_20210107_withBckAndValues.mat';
            anatomicalModel = load( sprintf('%s%s',basePath,fileName) );
        end
    case 8 % voxelman
         ME = MException('Model:wrongCourseID',...
        '%s : Voxelman model not available',functionName);
        throw(ME);    
    otherwise % MIDA as default
        fileName = 'MIDA_20210110.mat';
        anatomicalModel = load(sprintf('%s%s',basePath,fileName));
end

% correct to 6 properties
if anatomicalModel.numProperties < 6
    anatomicalModel.numProperties = 6;
    anatomicalModel.tissueValues(anatomicalModel.numTissues,6) = 0.0;
end
%% TODO
% add additional data
anatomicalModel.b0          = 1.0;
anatomicalModel.mu          = 1.0;
anatomicalModel.numVoxels   = size(anatomicalModel.spatial,1);
% diffusion: values in the ballpark based on T1 and T2
dummyDiffusion = 1e-9*( 0.1 ...
    + 0.2*anatomicalModel.tissueValues(:,1) ...
    + 2*anatomicalModel.tissueValues(:,2) );
anatomicalModel.tissueDiff = zeros(anatomicalModel.numTissues,3);
anatomicalModel.tissueDiff(:,1) = 0.9*dummyDiffusion;
anatomicalModel.tissueDiff(:,2) = 1.0*dummyDiffusion;
anatomicalModel.tissueDiff(:,3) = 0.8*dummyDiffusion;

%% report
tTotal = toc(tTotal);
if isempty(anatomicalModel)
    ME = MException('Model:emptyModel',...
        '%s : error loading -- empty anatomicalModel',functionName);
    throw(ME);
else
    fprintf(fid, '\n%s : %s model with %d voxels loaded',...
        functionName, anatomicalModel.name, anatomicalModel.numVoxels );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
end