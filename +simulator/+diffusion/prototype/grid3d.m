function r = grid3d(x, y, z)
%%    Generates a 3D grid
% _________________________________________________________________________
%
%       Generates a 3D cartesian grid or coordinates
%
% _________________________________________________________________________
%
%% INPUT
%   x           positions of the x coordinates
%
%
%% OPTIONAL INPUT
%   y           positions of the y coordinates
%   z           positions of the z coordinates
%
%
%% OUTPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if nargin == 1
    y = x;
    z = x;
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
