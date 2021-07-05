%%   Onboarding example 2
%    Bloch Simulation: Recieving signal integration
% 
% -------------------------------------------------------------------------
%
%    The goal of this script is to show and help to understand
%    the mechanics of the echo generation and implications in 
%    the signal received.
%
%    A uniform volume of 1mm3 is discretized in the X direction
%    and a varifying Gradient is applied.
%
%    Exercises and questions are posed at the end of the script to 
%    understand issues such as discretization and sputious echoes,
%    and how they related to situations in real life MRI. 
%
% -------------------------------------------------------------------------
%    REFERENCES:
%       https://en.wikipedia.org/wiki/Bloch_equations
%       https://drive.google.com/drive/folders/1xn2fsl06gcg5d0netH4yMKDfH6DFvWC9?usp=sharing
%
% -------------------------------------------------------------------------
%%  CORSMED AB © 2020
%   Jorge Fernandez Villena - jorge.villena@corsmed.com
% _________________________________________________________________________

clc;
clear all;
close all;

figure(1); clf;
title('ADJUST THE FIGURE SIZE TO YOUR SCREEN, THEN PRESS ENTER');
pause();

%% PARAMETERS TO PLAY WITH
NVOXELS = 5;        % number of voxels in the discretization
GSTRENGTH = 0.0e-3; % gradient strength
XSHIFT = 0.0e-3;    % shift of the voxel positions to create a non-symmetric problem
% Time constants
T1 = 0.800;  % longitudinal recovery time
T2 = 0.300;  % transversal decay
% initial magnetization: 
%    make sure magnitude is 1 to be consistent with MEQ (equilibrium)
Mox = 1.0;  % initaial x component of the magnetization
Moy = 0.0;  % initaial y component of the magnetization
Moz = 0.0;  % initaial z component of the magnetization
%% SET UP GRADIENT SIGNAL:
% generate gradient time signal:
DURATION = 200e-3;
TIMESTEP = 0.5e-3;
t = 0:1:ceil(DURATION/TIMESTEP);  % discretization of the time in integer steps
t = t*TIMESTEP;                   % time discretization in time steps 
% Gradient signal 
%   Zero for 10ms
%   Positive for 25ms
idxPlus1 = find( (10e-3 < t) & (t <= 25e-3) );
%   Zero for 15ms
%   Negative for (50ms) to generate a gradient echo at t=75ms
idxMinus1 = find( (50e-3 < t) & (t <= 100e-3) );
%   Positive for the remaining
idxPlus2 = find(t > 100e-3);
g = 0*t;
g(idxPlus1)  =  1;
g(idxMinus1) = -1;
g(idxPlus2)  =  1;

%% BASIC DEFINITIONS
BSTRENGTH = 1.5;        % main field strength
GAMMA = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
WGAMMA = 2*pi*GAMMA;    % angular gyromagnetic ratio

%% DOMAIN DEFINITIONS
FOV = 1e-3;
TOTALVOL = FOV.^3; % assume the total volume is FOV^3
DX  = FOV/NVOXELS;
X = linspace((-FOV+DX)/2, (FOV-DX)/2, NVOXELS).'; % spatial discretization
X = X + XSHIFT;   % shift the voxel positions to create a non-symmetric problem
VOL = DX*FOV*FOV; % with the discretization each voxels will have this volume

%% DEFINE TISSUE PROPERTIES
PD = 1.000*ones(NVOXELS,1);  % proton density: constant

%% INITIAL MAGNETIZATION
Mx = Mox*ones(NVOXELS,1);
My = Moy*ones(NVOXELS,1);
Mz = Moz*ones(NVOXELS,1);

MEQ = 1.0; % equilibrium magnetization in the Z direction

% transform T1 and T2 into vectors
T1 = T1*ones(NVOXELS,1);
T2 = T2*ones(NVOXELS,1);

% initalize variables for simulation and visualization
tprev = 0;   % initialize time
sx_time = 0*t;
sy_time = 0*t;
sz_time = 0*t;

for tt = 1:length(t)
    
    % compute current time step
    dt = t(tt) - tprev;
    
    % at each time step, assume constant fields: evaluate each component
    Bz = g(tt)*GSTRENGTH*X;  % z component at each voxel due to the gradient
  
    % Rotation due to Bz field
    zPhase = WGAMMA*Bz*dt; % phase accrued in the time step
    % apply rotation: note that these are now vectors (dot product)
    Mxtemp =  cos(zPhase).*Mx + sin(zPhase).*My;
    Mytemp = -sin(zPhase).*Mx + cos(zPhase).*My;
    % assign values
    Mx = Mxtemp;
    My = Mytemp;

    % transverse decay and longitudinal recovery to MEQ
    Mx = exp(-dt./T2).*Mx;
    My = exp(-dt./T2).*My;
    Mz = exp(-dt./T1).*Mz + MEQ*( 1 - exp(-dt./T1) );
    
    % update time
    tprev = t(tt);
    
    % apply signal integration (sum of the contributions of all voxels)
    sx_time(tt) = sum(VOL*PD.*Mx);
    sy_time(tt) = sum(VOL*PD.*My);
    sz_time(tt) = sum(VOL*PD.*Mz);
    
    % plot results
    figure(1);
    
    subplot(1,2,1); cla;
    hold on;
    % SCALE the signals by the total volume only for visual purposes
    plot(t(1:tt), g(1:tt), ': ', 'Color', '#D95319', 'LineWidth', 3); 
    plot(t(1:tt), sx_time(1:tt)/TOTALVOL, 'r-', 'LineWidth', 3); 
    plot(t(1:tt), sy_time(1:tt)/TOTALVOL, 'b-', 'LineWidth', 3);  
    plot(t(1:tt), sz_time(1:tt)/TOTALVOL, 'k-', 'LineWidth', 3); 
    xlim([0,DURATION]); ylim([-1.1*MEQ,1.1*MEQ]);
    legend('Gradient (normalized)', 'S_x', 'S_y', 'S_z');
    xlabel('Time (s)'); ylabel('Magnetization');
    title(sprintf('Signal components at t=%.3f ms', t(tt)*1e3));

    subplot(1,2,2);  cla;
    hold on;
    for vv = 1:NVOXELS
        plot3([0,Mx(vv)], [0,My(vv)], [0,Mz(vv)], 'Color', '#0072BD', 'LineWidth', 3);
        plot3(Mx(vv), My(vv), Mz(vv), 'o', 'MarkerFaceColor', '#D95319', 'LineWidth', 3);
    end
    xlim([-1.2*MEQ,1.2*MEQ]); ylim([-1.2*MEQ,1.2*MEQ]); zlim([-1.2*MEQ,1.2*MEQ]);
    xlabel('Mx'); ylabel('My'); zlabel('Mz'); grid on;
    title(sprintf('Voxels Magnetization vectors at t=%.3f ms', t(tt)*1e3));
    view(55,20);
    
    drawnow;
    
end

%% EXERCISE 1:
%%  Change the strength of the gradient GSTRENGTH to 0.5e-3
%%        and see what happens
%%  QUESTION: when the magnetization rephase, 
%%        is the signal peak larger, equal or smaller 
%%        than the signal with no gradient?
%%  HINT: Transverse Magnetization will dephase and rephase 
%%        with rotations at different speeds 
%%        (each voxel sees its own field)
%%
%% EXERCISE 2:
%%  Change the strength of the gradient GSTRENGTH to 1e-3 and 5e-3
%%        and see what happens
%%  QUESTION: do you think the received signals are accurate 
%%        as we increase the gradient strength?
%%  QUESTION: can a very large continous gradient be clinically useful?
%%  HINT: Dephasing will kill the transverse magnetization 
%%
%% EXERCISE 3:
%%  With increased GSTRENGTH, increase the number of voxels NVOXELS
%%        and see what happens.
%%  HINT: Better representation of the phase accross the domain will 
%%        remove some of the spurious echoes.
%%
%% EXERCISE 4:
%%  Shift the positions of the Voxels (XSHIFT = 0.5e-3) 
%%        and see what happens when you repeat the previous experiments.        
%%  HINT: Non-symmetric gradient distribution.
%%
