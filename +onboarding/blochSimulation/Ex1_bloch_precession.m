%%   Onboarding example 1
%    Bloch Simulation: Simulation and effect of Bz and precession
%
% -------------------------------------------------------------------------
%
%    The goal of this script is to show and help to understand
%    the effects of the magentic fields on a voxel magnetization.
%    It will show how the Bloch equations can be simulated if there
%    is no RF field.
%
%    A single voxel magnetization is simulated in the presence of 
%    a Z field. The effects of the field, as well as the T1 and T2
%    impact is shown.
%
%    Exercises and questions are posed at the end of the script to 
%    understand the impact of the different parameters on the 
%    evolution of the voxel magnetization.
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
GVALUE   = 0.0;  % constant value of the contribution to Bz field
% Time constants (purposedly low to shorten experiment)
T1 = 0.100;  % longitudinal recovery time
T2 = 0.030;  % transversal decay
% initial magnetization: 
%    make sure magnitude is 1 to be consistent with MEQ (equilibrium)
Mox = 1.0;  % initaial x component of the magnetization
Moy = 0.0;  % initaial y component of the magnetization
Moz = 0.0;  % initaial z component of the magnetization

%% BASIC DEFINITIONS
BSTRENGTH = 1.5;        % main field strength
GAMMA = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
WGAMMA = 2*pi*GAMMA;    % angular gyromagnetic ratio

% define the time step based on the precession frequency
% make sure that we have at least 10 samples per period
% but make sure is not too small and it takes forever to simulate
SPLRATE = 1/(20*GAMMA*GVALUE);
TIMESTEP = min(1e-3,SPLRATE);
TIMESTEP = max(TIMESTEP,1e-4);

%% DEFINE TISSUE PROPERTIES
PD = 1.000;  % proton density

%% SET UP RF SIGNAL: GAUSSIAN PROFILE
DURATION = 200e-3;   % duration of the sequence
duration = ceil(DURATION/TIMESTEP); % duration as integer multiple of timestep
t = 0:1:duration;  % discretization of the time in integer steps
t = t*TIMESTEP;       % time discretization in time steps 

%% INITIAL MAGNETIZATION
Mx = Mox;
My = Moy;
Mz = Moz;

MEQ = 1.0; % equilibrium magnetization in the Z direction

% initalize variables for simulation and visualization
tprev = 0;   % initialize time
mx_time = 0*t;
my_time = 0*t;
mz_time = 0*t;
bx_time = 0*t;
by_time = 0*t;
bz_time = 0*t;

for tt = 1:length(t)
    
    % compute current time step
    dt = t(tt) - tprev;
    
    % at each time step, assume constant fields: evaluate each component
    Bz = GVALUE;                  % z component of the field is gradient
  
    % Rotation due to Bz field
    zPhase = WGAMMA*Bz*dt; % phase accrued in the time step
    % generate rotation matrix
    A = [  cos(zPhase) , sin(zPhase) , 0 ; ...
          -sin(zPhase) , cos(zPhase) , 0 ; ...
           0           , 0           , 1 ];
    % apply rotation
    Mxtemp =  A(1,1)*Mx + A(1,2)*My + A(1,3)*Mz;
    Mytemp =  A(2,1)*Mx + A(2,2)*My + A(2,3)*Mz;
    Mztemp =  A(3,1)*Mx + A(3,2)*My + A(3,3)*Mz;
    % assign values
    Mx = Mxtemp;
    My = Mytemp;
    Mz = Mztemp;

    % transverse decay and longitudinal recovery to MEQ
    Mx = exp(-dt/T2)*Mx;
    My = exp(-dt/T2)*My;
    Mz = exp(-dt/T1)*Mz + MEQ*( 1 - exp(-dt/T1) );
    
    % update time
    tprev = t(tt);
    
    % plot results
    figure(1);
    
    % assign magnetizations to arrays
    mx_time(tt) = Mx;
    my_time(tt) = My;
    mz_time(tt) = Mz;

    subplot(1,2,1); cla;
    hold on;
    plot(t(1:tt), mx_time(1:tt), 'r-', 'LineWidth', 3); 
    plot(t(1:tt), my_time(1:tt), 'b-', 'LineWidth', 3);  
    plot(t(1:tt), mz_time(1:tt), 'k-', 'LineWidth', 3); 
    xlim([0,DURATION]); ylim([-1.2*MEQ,1.2*MEQ]);
    legend('Mx', 'My', 'Mz');
    xlabel('Time (s)'); ylabel('Magnetization');
    title(sprintf('Magnetization components at t=%.3f ms', t(tt)*1e3));

    subplot(1,2,2);  cla;
    plot3([0,Mx], [0,My], [0,Mz], 'Color', '#0072BD', 'LineWidth', 3);
    hold on;
    plot3(mx_time(1:tt), my_time(1:tt), mz_time(1:tt), '.', 'MarkerFaceColor', '#EDB120');
    plot3(Mx, My, Mz, 'o', 'MarkerFaceColor', '#D95319', 'LineWidth', 3); 
    xlim([-1.2*MEQ,1.2*MEQ]); ylim([-1.2*MEQ,1.2*MEQ]); zlim([-1.2*MEQ,1.2*MEQ]);
    xlabel('Mx'); ylabel('My'); zlabel('Mz'); grid on;
    title(sprintf('Magnetization vector at t=%.3f ms', t(tt)*1e3));
    view(55,20);
    
    drawnow;
    
end

%% EXERCISE 1:
%%  Change the time constants T1 and T2 and see the effect in the evolution of the magnetization
%%  HINT: T1 > T2 always, see the effect 
%%
%% EXERCISE 2:
%%  Change the GVALUE to constat small values (GVALUE=1e-6), 
%%        and see what happens as you increase the value (GVALUE=2e-6, GVALUE=5e-6)
%%  HINT: Z rotations are introduced due to Bz fields different of B0.
%%        This can be due to Gradients or B0 inhomogeneities.
%%        Gradients are using for spatial encodings.
%%
%% EXERCISE 3:
%%  Change the initial magnetization among Mx, My and Mz components 
%%        (make sure the magnitude remains 1, better done with rotations)
%%        and see what happens
%%  HINT: MRI sequences flip the Magnetization to certain position,  
%%        and then receive the XY components while applying a gradient.
%%        An interesting case is to put Mz=-1, Mx=My=0
%%