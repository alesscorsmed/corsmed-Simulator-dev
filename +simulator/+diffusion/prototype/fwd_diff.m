function [dMdx, dMdy, dMdz] = fwd_diff(M, nx, ny, nz, dx, dy, dz)
%%    Applies the forward differences
% _________________________________________________________________________
%
%       
% _________________________________________________________________________
%
%% INPUT
%
%
%% OUTPUT
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jorge.villena@corsmed.com
%   CORSMED AB
% _________________________________________________________________________

% scaling for spatial FD
scalex = 1/dx;
scaley = 1/dy;
scalez = 1/dz;

% compute forward differences
dMdx = 0*M;
dMdy = 0*M;
dMdz = 0*M;
for ii = 0:nx*ny*nz-1
    
    zid = floor(ii/(nx*ny));
    yid = floor((ii - zid*nx*ny)/(nx));
    xid = ii - zid*nx*ny - yid*nx;
    
    % correct for matlab indexes
    zid = zid + 1;
    yid = yid + 1;
    xid = xid + 1;
    idx = ii + 1;
    
    if (nx > 1)
        if (xid == 1)
            % fwd diff
            dMdx(idx) = scalex*(M(idx + 1) - M(idx));
        elseif (xid == nx)
            % bwd diff
            dMdx(idx) = scalex*(M(idx) - M(idx-1));
        else
            % default center diff
            dMdx(idx) = scalex*(M(idx + 1) - M(idx-1))/2;
        end
    end
    if (ny > 1)
        shift = nx;
        if (yid == 1)
            % fwd diff
            dMdy(idx) = scaley*(M(idx+shift) - M(idx));
        elseif (yid == ny)
            % bwd diff
            dMdy(idx) = scaley*(M(idx) - M(idx-shift));
        else
            % default center diff
            dMdy(idx) = scaley*(M(idx+shift) - M(idx-shift))/2;
        end
    end
    if (nz > 1)
        shift = nx*ny;
        if (zid == 1)
            % fwd diff
            dMdz(idx) = scalez*(M(idx+shift) - M(idx));
        elseif (zid == nz)
            % bwd diff
            dMdz(idx) = scalez*(M(idx) - M(idx-shift));
        else
            % default center diff
            dMdz(idx) = scalez*(M(idx+shift) - M(idx-shift))/2;
        end
    end
end
