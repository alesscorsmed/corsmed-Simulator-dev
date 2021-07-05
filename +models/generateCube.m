function [r,x,y,z] = generateCube(fovX,fovY,fovZ,dx,dy,dz)
%%    Generates a 3D grid for a Cube
%

% discretization
nx = ceil(fovX/dx);
ny = ceil(fovY/dy);
nz = ceil(fovZ/dz);
% Domain centered around 0
if (nx <= 1)
    x = 0;
else
    x = linspace(-fovX/2+dx/2,fovX/2-dx/2,nx);
end
if (ny <= 1)
    y = 0;
else
    y = linspace(-fovY/2+dy/2,fovY/2-dy/2,ny);
end
if (nz <= 1)
    z = 0;
else
    z = linspace(-fovZ/2+dz/2,fovZ/2-dz/2,nz);
end

% define the dimensions
L = length(x);
M = length(y);
N = length(z);

% allocate space
r = zeros(L,M,N,3);

% -------------------------------------------------------------------------
% Fill data
% -------------------------------------------------------------------------

for ix = 1:L
    xx = x(ix);
    for iy = 1:M
        yy = y(iy);
        for iz = 1:N
            zz = z(iz);
            r(ix,iy,iz,:) = [xx yy zz];
        end
    end
end
