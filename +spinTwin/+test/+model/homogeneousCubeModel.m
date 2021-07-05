function [spinModel] = homogeneousCubeModel(fovX,fovY,fovZ,dx,dy,dz,...
    t1Value, t2Value, b0, mu, debugMode, debugFile)
%
% MODELS.HOMOGENEOUSCUBEMODEL
%
%	Generates an homogenous cube model.
%
% INPUT
%   duration        total RF duration, in s
%   cycles          number of cycles of the sinc
%   angle           Flip angle, in degrees
%   tstep           time discretization
%   gamma           gyromagnetic ratio
%
% OUTPUT
%   time            discretized time vector, starts in tstep
%   signal          signal vector
%   BW              RF bandwidth, in Hz
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'models.homogeneousCubeModel';

if (nargin < 1 || isempty(fovX))
    fovX = 0.001;
end
if (nargin < 2 || isempty(fovY))
    fovY = 0.001; 
end
if (nargin < 3 || isempty(fovZ))
    fovZ = 0.001;
end
if (nargin < 4 || isempty(dx))
    dx = 1e-3;
end
if (nargin < 5 || isempty(dy))
    dy = 1e-3;
end
if (nargin < 6 || isempty(dz))
   dz = 1e-3;
end
if (nargin < 7 || isempty(t1Value))
    t1Value = 0.500;
end
if (nargin < 8 || isempty(t2Value))
    t2Value = 0.250;
end
if (nargin < 9 || isempty(b0))
    b0 = 1.5;
end
if (nargin < 10 || isempty(mu))
    mu = 1.0;
end
if (nargin < 11 || isempty(debugMode))
    debugMode = 0;
end
if (nargin < 12 || isempty(debugFile))
    debugFile = [];
end

% info for debugging
if debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(debugFile,'a');
    catch
        fid = 1;
    end
end

% generate a cube discretization
[r,x,y,z] = spinTwin.test.model.generateCube(fovX,fovY,fovZ,dx,dy,dz);
[nx,ny,nz,~] = size(r);
if nx > 1
    dx = x(2)-x(1);
else
    dx = fovX;
end
if ny > 1
    dy = y(2)-y(1);
else
    dy = fovY;
end
if nz > 1
    dz = z(2)-z(1);
else
    dz = fovZ;
end
niso = nx*ny*nz;

% % initialize spin model with 1 slice
% [spinModel] = data.spinModel.initialize(1);

% assign data
spinModel.name = 'Homogen Cube';
spinModel.numSlices = 1;
spinModel.totalIsochromats = niso;

spinModel.slice{1}.model.resolution     = [dx; dy; dz];
spinModel.slice{1}.model.numIsochromats = niso;
spinModel.slice{1}.model.nonZeroIndex   = 1:niso;

spinModel.slice{1}.model.r3D            = r;
spinModel.slice{1}.model.x              = reshape(r(:,:,:,1),nx*ny*nz,1);
spinModel.slice{1}.model.y              = reshape(r(:,:,:,2),nx*ny*nz,1);
spinModel.slice{1}.model.z              = reshape(r(:,:,:,3),nx*ny*nz,1);

spinModel.slice{1}.model.b0             = b0;
spinModel.slice{1}.model.mu             = mu;
spinModel.slice{1}.model.bi             = zeros(niso,1);
spinModel.slice{1}.model.pd             = ones(niso,1);


spinModel.slice{1}.model.numTissues         = 1; 
spinModel.slice{1}.model.numProperties      = 6; 
spinModel.slice{1}.model.tissueValues       = zeros(1,6);
spinModel.slice{1}.model.tissueValues(1,1)  = t1Value;
spinModel.slice{1}.model.tissueValues(1,2)  = t2Value;
spinModel.slice{1}.model.tissueType         = ones(niso,1);

numTissues = spinModel.slice{1}.model.numTissues;
spinModel.slice{1}.model.tissueDiff = zeros(numTissues,3);
%spinModel.slice{1}.model.xDiffusion = zeros(niso,1);
%spinModel.slice{1}.model.yDiffusion = zeros(niso,1);
%spinModel.slice{1}.model.zDiffusion = zeros(niso,1);

% coils
spinModel.slice{1}.model.numRxCoils  = 1;
spinModel.slice{1}.model.rxCoilMapsX = ones(niso,1);
spinModel.slice{1}.model.rxCoilMapsY = zeros(niso,1);

% final message
if  debugMode
    fprintf(fid, ...
        '\n%s done: %s generated with %d isochromats in %d slices\n',...
        functionName, spinModel.name, spinModel.totalIsochromats, spinModel.numSlices);
    if fid ~=1
        fclose(fid);
    end
end
