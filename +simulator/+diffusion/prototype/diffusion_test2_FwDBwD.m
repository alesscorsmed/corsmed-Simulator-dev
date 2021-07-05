%%    Pure diffusion test
%       Simple test, with homogeneous diffusion
%       Applies time domain simulation
%       Diffusion equation by simple Finite differences
%       Forward diff is applied for the gradient
%       Backward diff is applied for the div
%       Boundaries are set so that there is zero flux
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________

clear all;
close all;
%clc;
p = genpath(pwd);
addpath(p);

%% domain
% problem size and dimensions
fov_x = 0.100;
fov_y = 0.100;
fov_z = 0.001;
d_x = 1e-3;
d_y = 1e-3;
d_z = 1e-3;
n_x = fov_x/d_x;
n_y = fov_y/d_y;
n_z = fov_z/d_z;
% Domain centered around 0
if (n_x <= 1)
    x = 0;
else
    x = -n_x*d_x/2:d_x:(n_x/2)*d_x;
end
if (n_y <= 1)
    y = 0;
else
    y = -n_y*d_y/2:d_y:(n_y/2)*d_y;
end
if (n_z <= 1)
    z = 0;
else
    z = -n_z*d_z/2:d_z:(n_z/2)*d_z;
end
r = grid3d(x,y,z);
[n_x,n_y,n_z,~] = size(r);



% Diffussion terms: Positive direction defined
%   ASUMES that the term Dxx(i,j,k) is the diffusion from i,j,k -> i+1,j,k
%        diffusion i,j,k <- i+1,j,k is the same
%   ASUMES that the term Dyy(i,j,k) is the diffusion from i,j,k -> i,j+1,k
%   ASUMES that the term Dzz(i,j,k) is the diffusion from i,j,k -> i,j,k+1
% for the Cross-terms
%   Dxy(i,j,k) is cross diffusion from i,j,k -> i,j+1,k
%   Dyx(i,j,k) is cross diffusion from i,j,k -> i+1,j,k
% and so on
scale = 1e-6;
Dxx = scale*ones(n_x,n_y,n_z); % XX term
Dxy = scale*zeros(n_x,n_y,n_z);
Dxz = scale*zeros(n_x,n_y,n_z);
Dyx = scale*zeros(n_x,n_y,n_z);
Dyy = scale*ones(n_x,n_y,n_z); % YY term
Dyz = scale*zeros(n_x,n_y,n_z);
Dzx = scale*zeros(n_x,n_y,n_z);
Dzy = scale*zeros(n_x,n_y,n_z);
Dzz = scale*zeros(n_x,n_y,n_z);

% set up zeros for boundaries
Dxx(end,:,:) = 0;
Dyy(:,end,:) = 0;

% store spatial step for diffusion FD
scalex = 1/d_x;
scaley = 1/d_y;
scalez = 1/d_z;
vol = d_x*d_y*d_z;

% time steps
t = linspace(0,200,2001);

% initial distribution in xy plane
M = abs(r(:,:,:,1) + 1j*r(:,:,:,2));
M = (1 - M/max(M(:))).^3;

figure(1);
s = surf(M); s.EdgeColor = 'none'; % colormap hot;
zlim([0, 1]);
title(sprintf('t=%.3f, Sum=%g', t(1), vol*sum(M(:)) ));


fprintf(1, '\n initial volume %g \n', vol*sum(M(:)));

% time simulation
for tt = 2:length(t)
    
    % compute weighted Gradient by forward differences
    Vx = 0*M;
    Vy = 0*M;
    Vz = 0*M;
    for ii = 0:n_x*n_y*n_z-1
        
        zid = floor(ii/(n_x*n_y));
        yid = floor((ii - zid*n_x*n_y)/(n_x));
        xid = ii - zid*n_x*n_y - yid*n_x;
        
        % correct for matlab indexes
        zid = zid + 1;
        yid = yid + 1;
        xid = xid + 1;
        idx = ii + 1;
        
        if (n_x > 1) && (xid < n_x)
            % default Forward diff
            pi = idx + 1;
            dMdx = scalex*(M(pi) - M(idx));
        else % at end boundary
            dMdx = -scalex*M(idx);
        end
        if (n_y > 1) && (yid < n_y)
            % default Forward diff
            pi = idx + n_x;
            dMdy = scaley*(M(pi) - M(idx));
        else % at end boundary
            dMdy = -scaley*M(idx);
        end
        if (n_z > 1) && (zid < n_z)
            % default Forward diff
            pi = idx + n_x*n_y;
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
    for ii = 0:n_x*n_y*n_z-1
        
        zid = floor(ii/(n_x*n_y));
        yid = floor((ii - zid*n_x*n_y)/(n_x));
        xid = ii - zid*n_x*n_y - yid*n_x;
        
        % correct for matlab indexes
        zid = zid + 1;
        yid = yid + 1;
        xid = xid + 1;
        idx = ii + 1;
        
        if (n_x > 1) && (xid > 1 )
            % default Backward diff
            ni = idx - 1;
            dVdx = scalex*(Vx(idx) - Vx(ni));
        else % at start boundary
            dVdx = scalex*Vx(idx);
        end
        if (n_y > 1) && (yid > 1)
            % default Backward diff
            ni = idx - n_x;
            dVdy = scaley*(Vy(idx) - Vy(ni));
        else % at start boundary
            dVdy = scaley*Vy(idx);
        end
        if (n_z > 1) && (zid > 1)
            % default Backward diff
            ni = idx - n_x*n_y;
            dVdz = scalez*(Vz(idx) - Vz(ni));
        else % zero at start boundary
            dVdz = scalez*Vz(idx);
        end

        % add contribution
        dM(idx) = dVdx + dVdy + dVdz;
    end
    
    if tt == 2
        % Forward Euler time step
        h = t(tt) - t(tt-1);
        M = M + h*dM;
        dMm1 = dM;
    else % order 2
        h = t(tt) - t(tt-1);
        M = M + 3*h/2*dM - h/2*dMm1;
        dMm1 = dM;
    end
    
    figure(1);
    s = surf(M); s.EdgeColor = 'none';% colormap hot;
    zlim([0, 1]);
    title(sprintf('t=%.3f, Sum=%g', t(tt), vol*sum(M(:)) ));
end

fprintf(1, ' final volume %g \n', vol*sum(M(:)));
