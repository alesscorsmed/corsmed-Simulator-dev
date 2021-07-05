%%   Bloch-Torrey simulation
%
%    simple case of diffusion coupled with Bloch analytical
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________



clear all;
close all;
clc;
p = genpath('./diffusion');
addpath(p);

%% domain
% problem size and dimensions
fov_x = 0.100;  fov_y = 0.001;  fov_z = 0.001;
d_x = 10.0e-4;     d_y = 1e-3;     d_z = 1e-3;
% discretization
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


for Dscale = [1,3,5]*1e-9
    
    RESULTS{1} = [];
    testnum = 0;
    
    for Gstrength = [1,5,10]*0.033
        for Gscale = (0.2:0.2:1)
            for TAU = (20:10:100)*1e-3
  
                %% Sequence
                % encoding in x direction
                num_samples = 128; % artificially large to generate large gradient
                grad_amplitude = Gstrength; % scale to be narrow gradients
                grad_slewrate = 150*1000;
                tstep = 1e-6;
                [t_grad,s_grad,~] = numerical.sequences.generate_cartesian_encoding(fov_x,num_samples,grad_amplitude,grad_slewrate,tstep);
                
                %Dscale = 3e-9; % diffusion coefficient in m^2/s
                AG  = Gscale*Gstrength; % gradient amplitude
                TG  = t_grad(end); % gradient duration
                %TAU = 50e-3; % time between gradients (same as echo time)
                
                grad_bk = round(ceil( TG/1e-3 )*1e-3/tstep); % grad block, with pre and post time
                in_grad = round( TG/tstep ); % time of grad openend
                bf_grad = ceil( (grad_bk - in_grad)/2 ); % wait time before gradient
                af_grad = grad_bk - bf_grad - in_grad; % wait time after gradient
                
                tstart_fwgr = bf_grad;
                tstart_wait = grad_bk;
                tstart_rwgr = tstart_wait + round( TAU/tstep ) - grad_bk + bf_grad;
                tend = round( TAU/tstep ) + grad_bk + round( TAU/2/tstep );
                
                % global time
                t = 0:tend;
                t = t*tstep;
                % gradient signal
                g_signal = 0*t;
                g_signal((tstart_fwgr+1):(tstart_fwgr+length(s_grad))) =  AG*s_grad;
                g_signal((tstart_rwgr+1):(tstart_rwgr+length(s_grad))) = -AG*s_grad;
                % area:
                tdiff = [0, diff(t)];
                g_area = cumsum(tdiff.*g_signal);
                

                %% Diffussion terms: Positive direction defined
                %   ASUMES that the term Dxx(i,j,k) is the diffusion from i,j,k -> i+1,j,k
                %        diffusion i,j,k <- i+1,j,k is the same
                %   ASUMES that the term Dyy(i,j,k) is the diffusion from i,j,k -> i,j+1,k
                %   ASUMES that the term Dzz(i,j,k) is the diffusion from i,j,k -> i,j,k+1
                % for the Cross-terms
                %   Dxy(i,j,k) is cross diffusion from i,j,k -> i,j+1,k
                %   Dyx(i,j,k) is cross diffusion from i,j,k -> i+1,j,k
                % and so on
                Dxx = Dscale*ones(n_x,n_y,n_z); %(1 - (20*abs(r(:,:,:,1))).^2); %rand(n_x,n_y,n_z); % XX term
                Dxy = Dscale*zeros(n_x,n_y,n_z);
                Dxz = Dscale*zeros(n_x,n_y,n_z);
                Dyx = Dscale*zeros(n_x,n_y,n_z);
                Dyy = Dscale*zeros(n_x,n_y,n_z); % YY term
                Dyz = Dscale*zeros(n_x,n_y,n_z);
                Dzx = Dscale*zeros(n_x,n_y,n_z);
                Dzy = Dscale*zeros(n_x,n_y,n_z);
                Dzz = Dscale*zeros(n_x,n_y,n_z); % ZZ term
                
                % set up zeros for boundaries
                Dxx(end,:,:) = 0;
                Dyy(:,end,:) = 0;
                Dzz(:,:,end) = 0;
                
                % store spatial step for diffusion FD
                scalex = 1/d_x;
                scaley = 1/d_y;
                scalez = 1/d_z;
                vol = d_x*d_y*d_z;
                
                %% initialization
                
                mu = 1e-6; % magnetic susceptibility 1ppm
                GAMMA = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
                WGAMMA = 2*pi*GAMMA; % angular gyromagnetic ratio
                
                % tissue time constants
                t2 = 2.000;
                t1 = 3.000;
                R1 = 1/t1 * ones(n_x,n_y,n_z);
                R2 = 1/t2 * ones(n_x,n_y,n_z);
                
                % Gradient
                Gx = r(:,:,:,1);
                
                % magnetization definitions
                b0 = 1.5;
                Moz = b0 * ones(n_x,n_y,n_z);
                Mx = Moz; % assumes magnetization flipped to 90 at start
                My = 0*Mx;
                Mz = 0*Mx;
                
                % transverse
                Mm = abs(Mx + 1j*My);
                Mp = atan2(My,Mx);
                
                %% EXPONENTIAL DECAY DUE TO SELF DIFFUSION
                Bcoeff = (WGAMMA^2)*(AG^2)*(TG^2)*(TAU-TG/3); % s/m^2 1m2/1e6mm2
                ADC = exp(-Bcoeff*Dscale);
                
                figure(1000);
                plot(t,g_signal, 'LineWidth', 2);
                xlabel('time (s)'); ylabel('gradient (T)');
                title(sprintf('Signal with b=%.1g s/m^2, Gamp=%.2f mT, TE =%.2f ms Gt=%.2f ms ', Bcoeff, AG*1e3, TAU*1e3, TG*1e3));
                
                
                %% simulation
                Mx_df = Mx;
                My_df = My;
                Mz_df = Mz;
                
                % transverse
                Mm_df = abs(Mx_df + 1j*My_df);
                Mp_df = atan2(My_df,Mx_df);
                total_dMp = 0*Mp_df;
                
                
%                 % % create the video writer with 1 fps
%                 % writerObj = VideoWriter('diffusion_transverse_resolution_1mm.avi');
%                 % writerObj.FrameRate = 5;
%                 % open(writerObj);
%                 
%                 figure(10);
%                 clf(10);
%                 
%                 subplot(3,1,1);
%                 plot(x, Mm(:),'LineWidth', 2);
%                 hold on;
%                 plot(x, Mm_df(:),'LineWidth', 2);
%                 plot(x, ADC*Mm(:),'--', 'LineWidth', 1);
%                 legend('No Diffusion', 'w/ Diffusion', 'DW');
%                 xlabel('X position (m)'); ylabel('|M|');
%                 ylim([0 1.5])
%                 
%                 title( sprintf('Diffusion effect at %.1f mm resolution\n time %.3f ms\n(No Diffusion) Signal = %f\n(W/ Diffusion) Signal = %f\n(W/ ADC) Signal = %f',...
%                     d_x*1e3, t(1)*1e3, abs(Mx(51) + 1j*My(51)), abs(Mx_df(51) + 1j*My_df(51)), ADC*Mm(51) ) );
%                 
%                 subplot(3,1,2);
%                 plot(x, Mp(:),'LineWidth', 2);
%                 hold on;
%                 plot(x, Mp_df(:),'LineWidth', 2);
%                 legend('No Diffusion', 'w/ Diffusion');
%                 xlabel('X position (m)'); ylabel('Phase');
%                 
%                 subplot(3,1,3);
%                 plot(t,g_signal, 'b-', 'LineWidth', 2);
%                 hold on;
%                 stem(t(1),g_signal(1),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
%                 xlabel('time (s)'); ylabel('x gradient amplitude');
%                 
%                 % % convert the image to a frame
%                 % frame = getframe(gcf);
%                 % drawnow;
%                 % writeVideo(writerObj, frame);
%                 
%                 
%                 %
%                 % figure(1);
%                 % CLIM = [-2*b0, 2*b0];
%                 % subplot(2,2,1); imagesc(Mx.'); title('Mx, no diffusion'); colorbar;
%                 % subplot(2,2,2); imagesc(My.'); title('My, no diffusion'); colorbar;
%                 % subplot(2,2,3); imagesc(Mx_df.'); title('Mx, w/ diffusion'); colorbar;
%                 % subplot(2,2,4); imagesc(My_df.'); title('My, w/ diffusion'); colorbar;
%                 %
%                 % figure(2);
%                 % plot(t,g_signal, 'b-', 'LineWidth', 2);
%                 % hold on;
%                 % stem(t(1),g_signal(1),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
%                 % xlabel('time (s)'); ylabel('gradient amplitude');
%                 %
                
                % time simulation
                stride = 500;
                tm1 = 0;
                am1 = 0;
                for tt = 2:stride:length(t)-1
                    
                    % apply time step
                    h = t(tt) - tm1;
                    area = g_area(tt) - am1;
                    if h > 0
                        % accrued phase on the step
                        theta = WGAMMA*area*Gx;
                        
                        % Non-diffusion
                        % relaxation
                        Mm = exp(-R2*h).*Mm;
                        Mz = Moz + exp(-R1*h).*(Mz - Moz);
                        % precession
                        Mp = Mp + theta;
                        % components
                        Mx = Mm.*cos(Mp);
                        My = Mm.*sin(Mp);
                        
                        % Diffusion of magnitude and phase
                        dMm = diffusion(Mm_df, n_x, n_y, n_z, d_x, d_y, d_z,...
                            Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy);
                        dMp = diffusion(Mp_df, n_x, n_y, n_z, d_x, d_y, d_z,...
                            Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy);
                        % derivatives of magnitude and phase
                        [dMdx, dMdy, dMdz] = fwd_diff(Mm_df, n_x, n_y, n_z, d_x, d_y, d_z);
                        [dPdx, dPdy, dPdz] = fwd_diff(Mp_df, n_x, n_y, n_z, d_x, d_y, d_z);
                        % second order cross terms
                        d2dx = Dxx.*(-1*Mm_df.*dPdx.^2 - 1j*2*dPdx.*dMdx);
                        
%                         figure(20);
%                         clf(20);
%                         
%                         subplot(3,2,1);
%                         plot(x, dMm(:),'LineWidth', 2);
%                         xlabel('X position (m)'); ylabel('Mag diffusion');
%                         subplot(3,2,2);
%                         plot(x, dMp(:).*Mm_df(:),'LineWidth', 2);
%                         xlabel('X position (m)'); ylabel('Phase diffusion');
%                         subplot(3,2,3);
%                         plot(x, real(d2dx),'LineWidth', 2);
%                         xlabel('X position (m)'); ylabel('Real 2nd order');
%                         subplot(3,2,4);
%                         plot(x, imag(d2dx),'LineWidth', 2);
%                         xlabel('X position (m)'); ylabel('Imag 2nd order');
                        
                        Mreal = Mm_df + h*( real(d2dx) + dMm );
                        Mimag = h*( imag(d2dx) - dMp.*Mm_df );
                        dMm = Mm_df;
                        Mm_df = sqrt(Mreal.*Mreal + Mimag.*Mimag);
                        dMm = Mm_df - dMm;
                        dMp = atan2(Mimag,Mreal);
                        Mp_df = Mp_df + dMp;
                        
%                         subplot(3,2,5);
%                         plot(x, dMm,'LineWidth', 2);
%                         xlabel('X position (m)'); ylabel('Mag increment');
%                         subplot(3,2,6);
%                         plot(x, dMp,'LineWidth', 2);
%                         xlabel('X position (m)'); ylabel('Phase increment');
                        
                        
                        %         Mx_df = Mm_df.*cos(Mp_df + h*dMp) + h*dMm.*cos(Mp_df);
                        %         My_df = Mm_df.*sin(Mp_df + h*dMp) + h*dMm.*sin(Mp_df);
                        %         Mm_df = abs(Mx_df + 1j*My_df);
                        %         Mp_df = atan2(My_df,Mx_df);
                        % z component
                        dM = diffusion(Mz_df, n_x, n_y, n_z, d_x, d_y, d_z,...
                            Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy);
                        
                        Mz_df = Mz_df + h*dM;
                        
                        % relaxation with diffusion
                        Mm_df = exp(-R2*h).*Mm_df;
                        Mz_df = Moz + exp(-R1*h).*(Mz_df - Moz);
                        % precession with diffusion
                        Mp_df = Mp_df + theta;
                        % components
                        Mx_df = Mm_df.*cos(Mp_df);
                        My_df = Mm_df.*sin(Mp_df);
                        
                        % update previous time and area
                        tm1 = t(tt);
                        am1 = g_area(tt);
                    end
                    
                    
%                     figure(10);
%                     clf(10);
%                     
%                     subplot(3,1,1);
%                     plot(x, Mm(:),'LineWidth', 2);
%                     hold on;
%                     plot(x, Mm_df(:),'LineWidth', 2);
%                     plot(x, ADC*Mm(:),'--', 'LineWidth', 1);
%                     legend('No Diffusion', 'w/ Diffusion', 'DW');
%                     xlabel('X position (m)'); ylabel('|M|');
%                     ylim([0 1.5])
%                     
%                     title( sprintf('Diffusion effect at %.1f mm resolution\n time %.3f ms\n(No Diffusion) Signal = %f\n(W/ Diffusion) Signal = %f\n(W/ ADC) Signal = %f',...
%                         d_x*1e3, t(tt)*1e3, abs(Mx(51) + 1j*My(51)), abs(Mx_df(51) + 1j*My_df(51)), ADC*Mm(51) ) );
%                     
%                     subplot(3,1,2);
%                     plot(x, Mp(:),'LineWidth', 2);
%                     hold on;
%                     plot(x, Mp_df(:),'LineWidth', 2);
%                     legend('No Diffusion', 'w/ Diffusion');
%                     xlabel('X position (m)'); ylabel('Phase');
%                     
%                     subplot(3,1,3);
%                     plot(t,g_signal, 'b-', 'LineWidth', 2);
%                     hold on;
%                     stem(t(tt),g_signal(tt),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
%                     xlabel('time (s)'); ylabel('x gradient amplitude');
%                     %
%                     %     % convert the image to a frame
%                     %     frame = getframe(gcf);
%                     %     drawnow;
%                     %     writeVideo(writerObj, frame);
%                     
%                     %     figure(1);
%                     %     clf(1);
%                     %     CLIM = [-2*b0, 2*b0];
%                     %     subplot(2,2,1); imagesc(Mx.'); title('Mx, no diffusion'); colorbar;
%                     %     subplot(2,2,2); imagesc(My.'); title('My, no diffusion'); colorbar;
%                     %     subplot(2,2,3); imagesc(Mx_df.'); title('Mx, w/ diffusion'); colorbar;
%                     %     subplot(2,2,4); imagesc(My_df.'); title('My, w/ diffusion'); colorbar;
%                     %
%                     %     figure(2);
%                     %     clf(2);
%                     %     plot(t,g_signal, 'b-', 'LineWidth', 2);
%                     %     hold on;
%                     %     stem(t(tt),g_signal(tt),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
%                     %     xlabel('time (s)'); ylabel('gradient amplitude');
                    
                end
                
                % apply last time step
                h = t(end) - tm1;
                area = g_area(end) - am1;
                if h > 0
                    % accrued phase on the step
                    theta = WGAMMA*area*Gx;
                    
                    % Non-diffusion
                    % relaxation
                    Mm = exp(-R2*h).*Mm;
                    Mz = Moz + exp(-R1*h).*(Mz - Moz);
                    % precession
                    Mp = Mp + theta;
                    % components
                    Mx = Mm.*cos(Mp);
                    My = Mm.*sin(Mp);
                    
                    % Diffusion
                    % Diffusion of magnitude and phase
                    dMm = diffusion(Mm_df, n_x, n_y, n_z, d_x, d_y, d_z,...
                        Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy);
                    dMp = diffusion(Mp_df, n_x, n_y, n_z, d_x, d_y, d_z,...
                        Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy);
                    % derivatives of magnitude and phase
                    [dMdx, dMdy, dMdz] = fwd_diff(Mm_df, n_x, n_y, n_z, d_x, d_y, d_z);
                    [dPdx, dPdy, dPdz] = fwd_diff(Mp_df, n_x, n_y, n_z, d_x, d_y, d_z);
                    % second order cross terms
                    d2dx = Dxx.*(-1*Mm_df.*dPdx.^2 - 1j*2*dPdx.*dMdx);
                    
%                     figure(20);
%                     clf(20);
%                     
%                     subplot(3,2,1);
%                     plot(x, dMm(:),'LineWidth', 2);
%                     xlabel('X position (m)'); ylabel('Mag diffusion');
%                     subplot(3,2,2);
%                     plot(x, dMp(:).*Mm_df(:),'LineWidth', 2);
%                     xlabel('X position (m)'); ylabel('Phase diffusion');
%                     subplot(3,2,3);
%                     plot(x, real(d2dx),'LineWidth', 2);
%                     xlabel('X position (m)'); ylabel('Real 2nd order');
%                     subplot(3,2,4);
%                     plot(x, imag(d2dx),'LineWidth', 2);
%                     xlabel('X position (m)'); ylabel('Imag 2nd order');
                    
                    Mreal = Mm_df + h*( real(d2dx) + dMm );
                    Mimag = h*( imag(d2dx) - dMp.*Mm_df );
                    dMm = Mm_df;
                    Mm_df = sqrt(Mreal.*Mreal + Mimag.*Mimag);
                    dMm = Mm_df - dMm;
                    dMp = atan2(Mimag,Mreal);
                    Mp_df = Mp_df + dMp;
                    
%                     subplot(3,2,5);
%                     plot(x, dMm,'LineWidth', 2);
%                     xlabel('X position (m)'); ylabel('Mag increment');
%                     subplot(3,2,6);
%                     plot(x, dMp,'LineWidth', 2);
%                     xlabel('X position (m)'); ylabel('Phase increment');
                    
                    %         Mx_df = Mm_df.*cos(Mp_df + h*dMp) + h*dMm.*cos(Mp_df);
                    %         My_df = Mm_df.*sin(Mp_df + h*dMp) + h*dMm.*sin(Mp_df);
                    %         Mm_df = abs(Mx_df + 1j*My_df);
                    %         Mp_df = atan2(My_df,Mx_df);
                    % z component
                    dM = diffusion(Mz_df, n_x, n_y, n_z, d_x, d_y, d_z,...
                        Dxx, Dyy, Dzz, Dxy, Dyx, Dxz, Dzx, Dyz, Dzy);
                    
                    Mz_df = Mz_df + h*dM;
                    
                    % relaxation with diffusion
                    Mm_df = exp(-R2*h).*Mm_df;
                    Mz_df = Moz + exp(-R1*h).*(Mz_df - Moz);
                    % precession with diffusion
                    Mp_df = Mp_df + theta;
                    % components
                    Mx_df = Mm_df.*cos(Mp_df);
                    My_df = Mm_df.*sin(Mp_df);
                    
                end
                
%                 figure(10);
%                 clf(10);
%                 
%                 subplot(3,1,1);
%                 plot(x, Mm(:),'LineWidth', 2);
%                 hold on;
%                 plot(x, Mm_df(:),'LineWidth', 2);
%                 plot(x, ADC*Mm(:),'--', 'LineWidth', 1);
%                 legend('No Diffusion', 'w/ Diffusion', 'DW');
%                 xlabel('X position (m)'); ylabel('|M|');
%                 ylim([0 1.5])
%                 
%                 title( sprintf('Diffusion effect at %.1f mm resolution\n time %.3f ms\n(No Diffusion) Signal = %f\n(W/ Diffusion) Signal = %f\n(W/ ADC) Signal = %f',...
%                     d_x*1e3, t(end)*1e3, abs(Mx(51) + 1j*My(51)), abs(Mx_df(51) + 1j*My_df(51)), ADC*Mm(51) ) );
%                 
%                 subplot(3,1,2);
%                 plot(x, Mp(:),'LineWidth', 2);
%                 hold on;
%                 plot(x, Mp_df(:),'LineWidth', 2);
%                 legend('No Diffusion', 'w/ Diffusion');
%                 xlabel('X position (m)'); ylabel('Phase');
%                 
%                 subplot(3,1,3);
%                 plot(t,g_signal, 'b-', 'LineWidth', 2);
%                 hold on;
%                 stem(t(end),g_signal(end),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
%                 xlabel('time (s)'); ylabel('x gradient amplitude');
%                 
%                 % % convert the image to a frame
%                 % frame = getframe(gcf);
%                 % drawnow;
%                 % writeVideo(writerObj, frame);
%                 %
%                 %
%                 % % close the writer object
%                 % close(writerObj);
                
                
                
                figure(11);
                clf(11);
                
                subplot(3,1,1);
                plot(x, Mx(:),'LineWidth', 2);
                hold on;
                plot(x, Mx_df(:),'LineWidth', 2);
                plot(x, ADC*Mm(:),'--', 'LineWidth', 1);
                legend('No Diffusion', 'w/ Diffusion', 'DW');
                xlabel('X position (m)'); ylabel('Mx');
                
                title( sprintf('Diffusion effect at %.1f mm resolution\n time %.3f ms\n(No Diffusion) Signal = %f\n(W/ Diffusion) Signal = %f\n(W/ ADC) Signal = %f',...
                    d_x*1e3, t(end)*1e3, d_x*sum(abs(Mx(:) + 1j*My(:))), d_x*sum(abs(Mx_df(:) + 1j*My_df(:))),  d_x*sum(ADC*Mm(:)) ) );
                
                subplot(3,1,2);
                plot(x, My(:),'LineWidth', 2);
                hold on;
                plot(x, My_df(:),'LineWidth', 2);
                legend('No Diffusion', 'w/ Diffusion');
                xlabel('X position (m)'); ylabel('My');
                
                subplot(3,1,3);
                plot(t,g_signal, 'b-', 'LineWidth', 2);
                hold on;
                stem(t(end),g_signal(end),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
                xlabel('time (s)'); ylabel('x gradient amplitude');
                
                
                
                % norm(Mx_df_str1(:)) - norm(Mx_df(:))
                % norm(My_df_str1(:)) - norm(My_df(:))
                
                %
                % figure(1);
                % clf(1);
                % CLIM = [-2*b0, 2*b0];
                % subplot(2,2,1); imagesc(Mx.'); title('Mx, no diffusion'); colorbar;
                % subplot(2,2,2); imagesc(My.'); title('My, no diffusion'); colorbar;
                % subplot(2,2,3); imagesc(Mx_df.'); title('Mx, w/ diffusion'); colorbar;
                % subplot(2,2,4); imagesc(My_df.'); title('My, w/ diffusion'); colorbar;
                %
                % figure(2);
                % clf(2);
                % plot(t,g_signal, 'b-', 'LineWidth', 2);
                % hold on;
                % stem(t(end),g_signal(end),'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','red')
                % xlabel('time (s)'); ylabel('gradient amplitude');
                %
                %
                % figure(10);
                % clf(10);
                % CLIM = [-2*b0, 2*b0];
                % subplot(2,2,1); imagesc(Mx_df_str1.'); title('Mx, 1mm'); colorbar;
                % subplot(2,2,2); imagesc(My_df_str1.'); title('My, 1mm'); colorbar;
                % subplot(2,2,3); imagesc(Mx_df.'); title('Mx, 0.5mm'); colorbar;
                % subplot(2,2,4); imagesc(My_df.'); title('My, 0.5mm'); colorbar;
                %
                % figure(11);
                % clf(11);
                % CLIM = [-2*b0, 2*b0];
                % subplot(2,2,1); imagesc(abs(Mx_df_str1.' + 1j*My_df_str1.')); title('|M|, 1mm'); colorbar;
                % subplot(2,2,2); imagesc(angle(Mx_df_str1.' + 1j*My_df_str1.')); title('Phase, 1mm'); colorbar;
                % subplot(2,2,3); imagesc(abs(Mx_df.' + 1j*My_df.')); title('|M|, 0.5mm'); colorbar;
                % subplot(2,2,4); imagesc(angle(Mx_df.' + 1j*My_df.')); title('Phase, 0.5mm'); colorbar;
                %
                
                localResults.tstep = stride*tstep;
                localResults.resolution = d_x;
                localResults.Gs = Gstrength;
                localResults.Dc = Dscale;
                localResults.Ga = AG;
                localResults.Gt = TG;
                localResults.Te = TAU;
                localResults.Bc = Bcoeff;
                localResults.Sa = ADC;
                localResults.Ss = Mm_df(round(n_x/2))/Mm(round(n_x/2));
                testnum = testnum+1;
                RESULTS{testnum} = localResults;
                
                Barray(testnum) = Bcoeff;
            end
        end
    end
    
    figure(100+round(Dscale*1e9));
    bval = linspace(0,7,1001)*1e8;
    plot(bval, exp(-Dscale*bval), 'LineWidth', 2);
    hold on;
    for kk = 1:testnum
        localResults = RESULTS{kk};
        %plot(localResults.Bc, localResults.Sa, 's', 'LineWidth', 2, 'MarkerSize',10,'MarkerEdgeColor','#0072BD');
        plot(localResults.Bc, localResults.Ss, 'o', 'LineWidth', 2, 'MarkerSize',6,'MarkerEdgeColor','#D95319');
        plot(localResults.Bc, abs(localResults.Sa-localResults.Ss), 's', 'LineWidth', 2, 'MarkerSize',6,'MarkerEdgeColor','#A2142F');

        
%         switch localResults.Te
%             case 20*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, 'o', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, 'o', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 30*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, 's', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, 's', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 40*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, 'd', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, 'd', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 50*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, '^', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, '^', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 60*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, 'v', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, 'v', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 70*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, '<', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, '<', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 80*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, '>', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, '>', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             case 90*1e-3
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, 'h', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, 'h', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%             otherwise
%                 figure(100+round(Dscale*1e9));
%                 hold on;
%                 plot(localResults.Bc, localResults.Sa, 'p', 'MarkerSize',12,'MarkerEdgeColor','#0072BD');
%                 plot(localResults.Bc, localResults.Ss, 'p', 'MarkerSize',8,'MarkerEdgeColor','#D95319');
%         end
    end
    legend('Exp(-b*D)', 'Simulation', 'Error');
    xlabel('b (s/m^2)'); ylabel('Signal Decay');
    title(sprintf('Diffusion D=%.1g m^2/s\n Signal decay w.r.t. no Diffusion for different configurations',Dscale));
    grid on;

end
