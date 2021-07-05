%%   Onboarding example 3
%    Bloch Simulation: Simulation and effect of the RF
% 
% -------------------------------------------------------------------------
%
%    The goal of this script is to show and help to understand
%    the effects of the RF on a voxel magnetization.
%    It will show how the Bloch equations can be simulated when there 
%    is an RF field applied.
%    It will also show how the existence of a Z field during the RF 
%    affect the evolution of the voxel magnetization.
%
%    A single voxel magnetization is simulated in the presence of 
%    a RF and Z field. The effects of the fields is shown.
%
%    Exercises and questions are posed at the end of the script to 
%    understand the impact of the different parameters on the 
%    evolution of the voxel magnetization, including the area of the 
%    RF signal, the impact of the RF phase, the strength of the 
%    Z component, initial magnetization, ...
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
GVALUE   = 0.000;         % constant value of the contribution to Bz field
RF_AREA  = 5.8350e-09;    % area of the RF signal (pre-fixed for 90 degrees), directly related to the flip angle
RF_PHASE = 0.00;          % constant phase of the RF signal
% initial magnetization: 
%    make sure magnitude is 1 to be consistent with MEQ (equilibrium)
Mox = 0.0;  % initaial x component of the magnetization
Moy = 0.0;  % initaial y component of the magnetization
Moz = 1.0;  % initaial z component of the magnetization

%% BASIC DEFINITIONS
GAMMA = 42.577478518e6; % Hz⋅T−1, gyromagnetic ratio for 1H protons
WGAMMA = 2*pi*GAMMA;    % angular gyromagnetic ratio
TIMESTEP = 5e-6;        % define the time discretization step
BSTRENGTH = 1.5;        % main field strength

%% DEFINE TISSUE PROPERTIES
T1 = 0.500;  % longitudinal recovery time
T2 = 0.150;  % transversal decay
PD = 1.000;  % proton density

%% SET UP RF SIGNAL: GAUSSIAN PROFILE
RF_DURATION = 0.5e-3;   % duration of the sequence
RF_WIDTH = 0.075e-3;    % approximated "width" of the RF pulse
% compute the time signal
mu = RF_DURATION/2;
duration = ceil(RF_DURATION/TIMESTEP); % duration as integer multiple of timestep
t = 0:1:duration;  % discretization of the time in integer steps
t = t*TIMESTEP;       % time discretization in time steps 
s = RF_AREA/(RF_WIDTH*sqrt(2*pi))*exp(-(t - mu).^2 /(2*RF_WIDTH^2)); % signal time profile

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
    Bx = s(tt)*cos(RF_PHASE);     % apply the RF as x component of the transmit field
    By = s(tt)*sin(RF_PHASE);     % zero y component of the transmit field
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
    
    % rotation due to the Bxy fields
    Bxy = abs(Bx + 1j*By);   % effective field magnitude
    xyPhase = WGAMMA*Bxy*dt; % phase accrued in time step
    ax = Bx/Bxy; % auxiliar variable to compute the rotations
    ay = By/Bxy; % auxiliar variable to compute the rotations
    % generate rotation matrix
    A = [ (ax^2 + ay^2 *cos(xyPhase)) , (ax*ay)*(1-cos(xyPhase))    , -ay*sin(xyPhase) ; ...
          (ax*ay)*(1-cos(xyPhase))    , (ax^2 + ay^2 *cos(xyPhase)) ,  ax*sin(xyPhase) ; ...
           ay*sin(xyPhase)            , -ax*sin(xyPhase)            ,  cos(xyPhase)    ];
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
    
    bx_time(tt) = Bx;
    by_time(tt) = By;
    bz_time(tt) = Bz;
    
    subplot(1,3,1); cla;
    hold on;
    plot(t(1:tt), bx_time(1:tt)*1e6, 'r-', 'LineWidth', 3); 
    plot(t(1:tt), by_time(1:tt)*1e6, 'b-', 'LineWidth', 3);  
    plot(t(1:tt), bz_time(1:tt)*1e6, 'k-', 'LineWidth', 3);  
    legend('Bx', 'By', 'Bz');
    xlim([0,RF_DURATION]); ylim([-50,50]);
    xlabel('Time (s)'); ylabel('Field Strength (uT)');
    title(sprintf('Sequence Signals at t=%.3f ms', t(tt)*1e3));

    subplot(1,3,2); cla;
    hold on;
    plot(t(1:tt), mx_time(1:tt), 'r-', 'LineWidth', 3); 
    plot(t(1:tt), my_time(1:tt), 'b-', 'LineWidth', 3);  
    plot(t(1:tt), mz_time(1:tt), 'k-', 'LineWidth', 3); 
    xlim([0,RF_DURATION]); ylim([-1.2*MEQ,1.2*MEQ]);
    legend('Mx', 'My', 'Mz');
    xlabel('Time (s)'); ylabel('Magnetization');
    title(sprintf('Magnetization components at t=%.3f ms', t(tt)*1e3));

    subplot(1,3,3);  cla;
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
%%  Change the Area of the RF and see the different flip angles achieved
%%  HINT: flip angle should be close to linear with the RF area

%% EXERCISE 2:
%%  Change the RF phase to constant PI/2, PI and -PI/2 ( and any other ) and see the effect
%%  HINT: depending on how we distribute the RF on the x and y components, we flip around a different axis

%% EXERCISE 3:
%%  Change the GVALUE to constat small values (GVALUE=10e-6), 
%%        and see what happens as you increase the value (GVALUE=25e-6, GVALUE=50e-6)
%%  HINT: Z rotations are introduced due to Bz fields different of B0.
%%        This can be due to Gradients, B0 inhomogeneities or Off-resonance RF signals
%%        The effect is different flip angles, and different final Mxy phases

%% EXERCISE 4:
%%  Change the initial magnetization among Mx, My and Mz components (make sure the magnitude remains 1)
%%        and see what happens
%%  HINT: An important case is to see what happens with Mangetization in the XY plane
%%        and a 180 degree Flip Angle
%%        This is the basis for the Spin Echo based sequences
