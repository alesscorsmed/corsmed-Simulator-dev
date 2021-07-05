function [sbrModel] = generateSBRmodel(fovX, fovY, fovZ, nX, nY, nZ,...
    pdValue, t1Value, t2Value, b0, mu)
%
% 
%
%========================  CORSMED AB © 2020 ==============================
%

% discretization
dx = fovX/nX;
dy = fovY/nY;
dz = fovZ/nZ;
% Domain centered around 0
if (nX <= 1)
    x = 0;
else
    x = linspace(-fovX/2+dx/2,fovX/2-dx/2,nX);
end
if (nY <= 1)
    y = 0;
else
    y = linspace(-fovY/2+dy/2,fovY/2-dy/2,nY);
end
if (nZ <= 1)
    z = 0;
else
    z = linspace(-fovZ/2+dz/2,fovZ/2-dz/2,nZ);
end
% allocate space
r = zeros(nX,nY,nZ,3);
% fill positions
for ix = 1:nX
    xx = x(ix);
    for iy = 1:nY
        yy = y(iy);
        for iz = 1:nZ
            zz = z(iz);
            r(ix,iy,iz,:) = [xx yy zz];
        end
    end
end

% assign data
niso = nX*nY*nZ;
sbrModel.resolution     = [dx; dy; dz];
sbrModel.numIsochromats = niso;
sbrModel.nonZeroIndex   = 1:niso;
% positions
sbrModel.r3D            = r;
sbrModel.x              = reshape(r(:,:,:,1),nX*nY*nZ,1);
sbrModel.y              = reshape(r(:,:,:,2),nX*nY*nZ,1);
sbrModel.z              = reshape(r(:,:,:,3),nX*nY*nZ,1);
%data
sbrModel.b0             = b0;
sbrModel.mu             = mu;
sbrModel.bi             = zeros(niso,1);
sbrModel.cs             = zeros(niso,1);
sbrModel.pr             = pdValue*ones(niso,1);
sbrModel.pi             = zeros(niso,1);
sbrModel.r1             = (1/max(t1Value,1e-2))*ones(niso,1);
sbrModel.r2             = (1/max(t2Value,1e-3))*ones(niso,1);
% coils
sbrModel.numRxCoils  = 1;
sbrModel.rxCoilMapsX = ones(niso,1);
sbrModel.rxCoilMapsY = zeros(niso,1);

