function [dM] = diffusion(M, nx, ny, nz, dx, dy, dz,...
                          Dxx, Dyy, Dzz, ...
                          Dxy, Dyx, Dxz, Dzx, Dyz, Dzy)
%%    Applies the diffusion
% _________________________________________________________________________
%
%       Diffusion equation by simple Finite differences
%       Forward diff is applied for the gradient
%       Backward diff is applied for the div
%       Boundaries are set so that there is zero flux
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

% compute weighted Gradient by forward differences
Vx = 0*M;
Vy = 0*M;
Vz = 0*M;
for ii = 0:nx*ny*nz-1
    
    zid = floor(ii/(nx*ny));
    yid = floor((ii - zid*nx*ny)/(nx));
    xid = ii - zid*nx*ny - yid*nx;
    
    % correct for matlab indexes
    zid = zid + 1;
    yid = yid + 1;
    xid = xid + 1;
    idx = ii + 1;
    
    if (nx > 1) && (xid < nx)
        % default Forward diff
        pi = idx + 1;
        dMdx = scalex*(M(pi) - M(idx));
    else % at end boundary
        dMdx = -scalex*M(idx);
    end
    if (ny > 1) && (yid < ny)
        % default Forward diff
        pi = idx + nx;
        dMdy = scaley*(M(pi) - M(idx));
    else % at end boundary
        dMdy = -scaley*M(idx);
    end
    if (nz > 1) && (zid < nz)
        % default Forward diff
        pi = idx + nx*ny;
        dMdz = scalez*(M(pi) - M(idx));
    else
        dMdz = -scalez*M(idx);
    end
    
    % compute v: weighted gradient
    Vx(idx) = Dxx(idx).*dMdx + Dxy(idx).*dMdy + Dxz(idx).*dMdz;
    Vy(idx) = Dyx(idx).*dMdx + Dyy(idx).*dMdy + Dyz(idx).*dMdz;
    Vz(idx) = Dzx(idx).*dMdx + Dzy(idx).*dMdy + Dzz(idx).*dMdz;
    
end

% compute Divergence by backward differences
dM = 0*M;
for ii = 0:nx*ny*nz-1
    
    % corresponding x, y and z indexes, assuming x,y,z ordering
    zid = floor(ii/(nx*ny));
    yid = floor((ii - zid*nx*ny)/(nx));
    xid = ii - zid*nx*ny - yid*nx;
    
    % correct for matlab indexes
    zid = zid + 1;
    yid = yid + 1;
    xid = xid + 1;
    idx = ii + 1;
    
    if (nx > 1) && (xid > 1 )
        % default Backward diff
        ni = idx - 1;
        dVdx = scalex*(Vx(idx) - Vx(ni));
    else % at start boundary
        dVdx = scalex*Vx(idx);
    end
    if (ny > 1) && (yid > 1)
        % default Backward diff
        ni = idx - nx;
        dVdy = scaley*(Vy(idx) - Vy(ni));
    else % at start boundary
        dVdy = scaley*Vy(idx);
    end
    if (nz > 1) && (zid > 1)
        % default Backward diff
        ni = idx - nx*ny;
        dVdz = scalez*(Vz(idx) - Vz(ni));
    else % zero at start boundary
        dVdz = scalez*Vz(idx);
    end
    
    % add contribution
    dM(idx) = dVdx + dVdy + dVdz;
end
