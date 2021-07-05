function [DR, DI] = diffusion_pm(M, P, nx, ny, nz, dx, dy, dz,...
                                 Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy)
%%    Applies the diffusion based on phase and magnitude
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
phase_limit = pi/8;
grad_step = dx;
scalex = 1/grad_step;


% compute weighted Gradient by forward differences
DR = 0*M;
DI = 0*M;
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
        if (xid < nx)
            % default Forward diff
            % force spatial step s.t. delta phase = scalep
            % and compute the corresponding plus on magnitude
            idx_pls = idx + 1;
            spatial_scale = grad_step/dx;
            phase_pls = (P(idx_pls) - P(idx))*spatial_scale; % delta phase
            if abs(phase_pls) > phase_limit
                spatial_scale = spatial_scale*phase_limit/abs(phase_pls); % scaled spatial delta 
            end
            dx_pls = dx*spatial_scale;
            phase_pls = (P(idx_pls) - P(idx))*spatial_scale; % actual delta phase
            mag_pls = M(idx) + (M(idx_pls) - M(idx))*spatial_scale; % plus one mag
            scaled_diff_pls = Dxx(idx)/dx_pls;
        else
            mag_pls = 0;
            phase_pls = 0;
            scaled_diff_pls = 0;
            dx_pls = dx_mns;
        end
        if (xid > 1)
            % default backward diff
            % force spatial step s.t. delta phase = scalep
            % and compute the corresponding plus on magnitude
            idx_mns = idx - 1;
            spatial_scale = grad_step/dx;
            phase_mns = (P(idx_mns) - P(idx))*spatial_scale; % delta phase
            if abs(phase_mns) > phase_limit
                spatial_scale = spatial_scale*phase_limit/abs(phase_mns); % scaled spatial delta 
            end
            dx_mns = dx*spatial_scale;
            phase_mns = (P(idx_mns) - P(idx))*spatial_scale; % actual delta phase
            mag_mns = M(idx) + (M(idx_mns) - M(idx))*spatial_scale; % minus one mag
            scaled_diff_mns = Dxx(idx_mns)/dx_mns;
        else
            mag_mns = 0;
            phase_mns = 0;
            scaled_diff_mns = 0;
            dx_mns = dx_pls;
        end
        %scalex = 2/(dx_pls+dx_mns);
        % compute Diffusion
        DR(idx) = scalex*( scaled_diff_pls*( mag_pls*cos(phase_pls) - M(idx) )...
                         + scaled_diff_mns*( mag_mns*cos(phase_mns) - M(idx) ) );
        DI(idx) = scalex*( scaled_diff_mns*mag_mns*sin(phase_mns) ...
                         - scaled_diff_pls*mag_pls*sin(phase_pls) );
        if xid == 51
            fprintf(' dx+ = %.3g, dx- = %.3g,  DR = %.3g,  DI = %.3g, \n', dx_pls, dx_mns, dx*DR(idx), dx*DI(idx))
        end
    end   
end

figure(100);
clf(100);
plot(DR(:),'LineWidth', 2);
hold on;
plot(DI(:),'LineWidth', 2);
legend('Real', 'Imag');
xlabel('X position (m)'); ylabel('D component');



