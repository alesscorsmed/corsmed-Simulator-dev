function [dMdxx, dMdxy, dMdxz, dMdyx, dMdyy, dMdyz, dMdzx, dMdzy, dMdzz ]...
         = dweighted_differences(M, nx, ny, nz, dx, dy, dz,...
                                Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy)
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
dMdxx = 0*M; dMdxy = 0*M; dMdxz = 0*M;
dMdyx = 0*M; dMdyy = 0*M; dMdyz = 0*M;
dMdzx = 0*M; dMdzy = 0*M; dMdzz = 0*M;
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
        shift = 1;
        if (xid == 1)
            % fwd diff
            dMdxx(idx) = scalex*(Dxx(idx+shift)*M(idx+shift) - Dxx(idx)*M(idx));
            dMdyx(idx) = scalex*(Dyx(idx+shift)*M(idx+shift) - Dyx(idx)*M(idx));
            dMdzx(idx) = scalex*(Dzx(idx+shift)*M(idx+shift) - Dzx(idx)*M(idx));
        elseif (xid == nx)
            % bwd diff
            dMdxx(idx) = scalex*(Dxx(idx)*M(idx) - Dxx(idx-shift)*M(idx-shift));
            dMdyx(idx) = scalex*(Dyx(idx)*M(idx) - Dyx(idx-shift)*M(idx-shift));
            dMdzx(idx) = scalex*(Dzx(idx)*M(idx) - Dzx(idx-shift)*M(idx-shift));
        else
            % default center diff
            dMdxx(idx) = scalex*(Dxx(idx+shift)*M(idx+shift) - Dxx(idx-shift)*M(idx-shift))/2;
            dMdyx(idx) = scalex*(Dyx(idx+shift)*M(idx+shift) - Dyx(idx-shift)*M(idx-shift))/2;
            dMdzx(idx) = scalex*(Dzx(idx+shift)*M(idx+shift) - Dzx(idx-shift)*M(idx-shift))/2;
        end
    end
    if (ny > 1)
        shift = nx;
        if (yid == 1)
            % fwd diff
            dMdxy(idx) = scaley*(Dxy(idx+shift)*M(idx+shift) - Dxy(idx)*M(idx));
            dMdyy(idx) = scaley*(Dyy(idx+shift)*M(idx+shift) - Dyy(idx)*M(idx));
            dMdzy(idx) = scaley*(Dzy(idx+shift)*M(idx+shift) - Dzy(idx)*M(idx));
        elseif (yid == ny)
            % bwd diff
            dMdxy(idx) = scaley*(Dxy(idx)*M(idx) - Dxy(idx-shift)*M(idx-shift));
            dMdyy(idx) = scaley*(Dyy(idx)*M(idx) - Dyy(idx-shift)*M(idx-shift));
            dMdzy(idx) = scaley*(Dzy(idx)*M(idx) - Dzy(idx-shift)*M(idx-shift));
        else
            % default center diff
            dMdxy(idx) = scaley*(Dxy(idx+shift)*M(idx+shift) - Dxy(idx-shift)*M(idx-shift))/2;
            dMdyy(idx) = scaley*(Dyy(idx+shift)*M(idx+shift) - Dyy(idx-shift)*M(idx-shift))/2;
            dMdzy(idx) = scaley*(Dzy(idx+shift)*M(idx+shift) - Dzy(idx-shift)*M(idx-shift))/2;
        end
    end
    if (nz > 1)
        shift = nx*ny;
        if (zid == 1)
            % fwd diff
            dMdxz(idx) = scalez*(Dxz(idx+shift)*M(idx+shift) - Dxz(idx)*M(idx));
            dMdyz(idx) = scalez*(Dyz(idx+shift)*M(idx+shift) - Dyz(idx)*M(idx));
            dMdzz(idx) = scalez*(Dzz(idx+shift)*M(idx+shift) - Dzz(idx)*M(idx));
        elseif (zid == nz)
            % bwd diff
            dMdxz(idx) = scalez*(Dxz(idx)*M(idx) - Dxz(idx-shift)*M(idx-shift));
            dMdyz(idx) = scalez*(Dyz(idx)*M(idx) - Dyz(idx-shift)*M(idx-shift));
            dMdzz(idx) = scalez*(Dzz(idx)*M(idx) - Dzz(idx-shift)*M(idx-shift));
        else
            % default center diff
            dMdxz(idx) = scalez*(Dxz(idx+shift)*M(idx+shift) - Dxz(idx-shift)*M(idx-shift))/2;
            dMdyz(idx) = scalez*(Dyz(idx+shift)*M(idx+shift) - Dyz(idx-shift)*M(idx-shift))/2;
            dMdzz(idx) = scalez*(Dzz(idx+shift)*M(idx+shift) - Dzz(idx-shift)*M(idx-shift))/2;
        end
    end
end
