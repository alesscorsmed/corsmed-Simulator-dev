%% INITIAL CONFIGURATION
% Model-specific parameters
model_desired_sizeofelement = [0.001,0.001,0.0005];
model_desired_dim = [0.3, 0.2, 0.3];

gamma           = 42.56*10^6;   % Hz/T
dt              = 5*10^-6;        % sec
sliceThickness  = 0.006;        % desired slice thickness (in m)
f0              = 64000000;     % in Hz

fast            = 0;
save            = 0;
plotON          = 1;
a2prepulse      = 'ramp';

RFmatrix        = zeros(6,1);
RFmatrix(1,:)   = 0.000480;         % RF duration 
RFmatrix(2,:)   = 2;                % RF sinc cycles
RFmatrix(3,:)   = [60];             % RF rotation angle
RFmatrix(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
RFmatrix(5,:)   = [1];              % starting timestep of RF
RFmatrix(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis

receiverBW      = 100*10^3;  % in Hz (this is a feature of the hardware)
G_max           = 40*10^-3; % in T/m

kspace(1,1) = 128; %defines the size of k-space along x direction  %128x128
kspace(1,2) = 128; %defines the size of k-space along y direction

FOV(1,1) = 0.2;
FOV(1,2) = 0.2;

acquiredkspace      = kspace;
acquiredFOV         = FOV;

partialFourierStruct.type   = 'readConjugate';
partialFourierStruct.factor = 1;%0.625;

lengthRamp  = 10;
addTailAby2 = 1;

TE = 0.0014;
TR = 0.0028;        % sec

receiver_BW_star = 1/dt; % this is the BW which depends on the time step 
                        % which is chosen for the simulation. It is needed
                        % to calculate the factor_BW by which the kspace(1,1)
                        % will be multiplied. This product will define the
                        % points that will be collected from the signal
                        % during the acquisition.

[pulse_sequence_bSSFP,isInKspace_bSSFP,pulse_sequence_TR_bSSFP,...
    soft_crushers_bSSFP,N_pulse_bSSFP,N_TR_bSSFP,kspace_times_bSSFP,BW_rf,...
    TR_times_bSSFP,receiver_phase_bSSFP,TRmin_bSSFP,PulseqDuration_bSSFP,...
    bSSFPCenterOfAcquisition,encodedkSpace] = ...
    pulseSequenceGenerator.bSSFPgenerator(RFmatrix,...
    kspace,FOV,receiverBW,G_max,TR,TE,dt,gamma,a2prepulse,[],[],[],...
    partialFourierStruct,lengthRamp,addTailAby2);
  
%% ALARM - CHECK IF SPURIOUS ECHOES WILL APPEAR
% PHASE DIFFERENCE OF CONTIGUOUS ELEMENTS ALARM
elementsize=model_desired_sizeofelement;
alarm=[1,1,1];
alarm_count = 0;
while alarm(1,1)==1||alarm(1,2)==1||alarm(1,3)==1
    [elementsize,alarm] = ...
        pulseSequenceGenerator.testAccumPhaseGRE(N_pulse_bSSFP,...
        N_TR_bSSFP,dt,pulse_sequence_bSSFP(4:6,:),elementsize);
    alarm_count = alarm_count + 1;
end

disp(['The alarm was enabled ',num2str(alarm_count - 1),' times.'])
disp(['The elementsize should be ',num2str(elementsize(1,1)),' x ',...
    num2str(elementsize(1,2)),' x ',num2str(elementsize(1,3)),' meters (m)'])

%% FIX THE PULSE SEQUENCE TO WORK PROPERLY WITH coreMRI
info_bSSFP = pulseSequenceGenerator.createInfoStruct_bSSFP(kspace,FOV,dt,...
    receiverBW,TR,TE,kspace_times_bSSFP,f0,sliceThickness,...
    receiver_phase_bSSFP,RFmatrix,acquiredkspace,acquiredFOV,encodedkSpace);

pulse_sequence_READOUT.pulse_sequence   = pulse_sequence_bSSFP;
pulse_sequence_READOUT.isInKspace       = isInKspace_bSSFP;
pulse_sequence_READOUT.soft_crushers    = soft_crushers_bSSFP;
pulse_sequence_READOUT.N_pulse          = N_pulse_bSSFP;
pulse_sequence_READOUT.kspace_times     = kspace_times_bSSFP;
pulse_sequence_READOUT.TR_times         = TR_times_bSSFP;
pulse_sequence_READOUT.receiver_phase   = receiver_phase_bSSFP;

% Formulate the MOLLI pulse sequence
MOLLI_scheme            = [5,3];
pause_cc                = 3;  % number of heartbeats
cardiac_cycle_duration  = 1;

% TI_initial_1 holds the TriggerTime of the DICOM images of these MOLLIs
TI_initial              = [250, 1250, 2250, 3250, 4250, 290, 1290, 2290];

TD_initial              = 0.5;  % choose a low TD_initial to fit the bSSFP readout 
                                % within the cardiac cycle

% TD is the time from the R-wave triggering till the start of the ramp
TD                      = abs(repmat(TD_initial,1,size(TI_initial,2)));
dicomMaxRR              = 1.5*cardiac_cycle_duration;

% Inversion pulse
% [B1,Dph,Dfr] = pulseSequenceGenerator.addRF_ADIABATIC_TAN_TANH_FULL_PASSAGE(0.000015,...
%     9500,10,1.525373047373320,IR_duration,dt);
[B1,Dph,Dfr] = pulseSequenceGenerator.addRF_AdiabaticHypSec(0.0000270606,...
    1343.39,0.00474,0.000005);

pulse_sequence_IR.pulse_sequence   = [B1;Dph;Dfr;zeros(5,size(B1,2))];
pulse_sequence_IR.isInKspace       = zeros(1,size(B1,2));
pulse_sequence_IR.soft_crushers    = zeros(1,size(B1,2));
pulse_sequence_IR.N_pulse          = size(B1,2);
pulse_sequence_IR.kspace_times     = [];
pulse_sequence_IR.TR_times         = [];
pulse_sequence_IR.receiver_phase   = [];

% Create MOLLI pulse sequence
[pulse_sequence_MOLLI,times_fitting_single_point] = ...
    pulseSequenceGenerator.MOLLIgenerator_ImprovedPerformance(...
    pulse_sequence_READOUT,info_bSSFP,pulse_sequence_IR,...
    cardiac_cycle_duration,TD,TI_initial,dt,MOLLI_scheme,pause_cc,TE,TR,...
    RFmatrix(1,:),acquiredkspace,[],[],[]);

pulse_sequence  = pulse_sequence_MOLLI.pulse_sequence;
N_pulse         = pulse_sequence_MOLLI.N_pulse;
infoMOLLI       = pulse_sequence_MOLLI.info;
soft_crushers   = pulse_sequence_MOLLI.soft_crushers;
isInKspace      = pulse_sequence_MOLLI.isInKspace;
receiver_phase  = pulse_sequence_MOLLI.receiver_phase;
kspace_times    = pulse_sequence_MOLLI.kspace_times;

info = pulseSequenceGenerator.createInfoStruct_MOLLI(kspace,FOV,dt,receiverBW,TR,TE,...
    kspace_times,f0,sliceThickness,receiver_phase,RFmatrix,acquiredkspace,...
            acquiredFOV,encodedkSpace,infoMOLLI,TI_initial);

%% ADD FAST ALGORITHM
% if fast == 1
%     [pulse_sequence_new,kspace_times_new,isInKspace_new,soft_crushers_new] = ...
%         pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
%         info.pulseSequence.kspace,soft_crushers);
%     
%     info_new = pulseSequenceGenerator.createInfoStruct_bSSFP(kspace,FOV,dt,...
%         receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,receiver_phase,...
%         RFmatrix,acquiredkspace,acquiredFOV,encodedkSpace);
% end

%% SAVE PULSE SEQUENCE
if save == 1
    
    % Save also the TR_times in the info structure
    info.pulseSequence.TR_times = TR_times;
    
    pulseSequenceGenerator.savePulseSequence(pulse_sequence,...
        soft_crushers,isInKspace,N_pulse,info,gamma,dt,'bSSFP_128x128')
    
    if fast == 1
        pulseSequenceGenerator.savePulseSequence(pulse_sequence_new,...
            soft_crushers_new,isInKspace_new,size(pulse_sequence_new,2),...
            info_new,gamma,dt,'bSSFP_128x128_fast')
    end
    
end

%% PLOT
if plotON == 1
    timeaxis = [1:N_pulse]*dt;
    soft_crushers_ind = find(soft_crushers);

    a1 = subplot(7,1,1);
    plot(timeaxis,pulse_sequence(1,:),'r')
    hold on
    plot(timeaxis(soft_crushers_ind),pulse_sequence(1,soft_crushers_ind),'b*')
    ylabel('RF magnitude') % in Tesla

    a2 = subplot(7,1,2);
    plot(timeaxis,pulse_sequence(2,:),'r')
    ylabel('RF phase') % in radians

    a3 = subplot(7,1,3);
    plot(timeaxis,pulse_sequence(3,:),'r')
    ylabel('RF Freq.') % in Hz

    a4 = subplot(7,1,4);
    plot(timeaxis,pulse_sequence(4,:))
    hold on
    plot(timeaxis(1,kspace_times(:)),pulse_sequence(4,kspace_times(:)),'go')
    hold on
%     plot(timeaxis(1,times_fitting_single_point),pulse_sequence(4,times_fitting_single_point),'b*')
    ylabel('Gx') % in T/m

    a5 = subplot(7,1,5);
    plot(timeaxis,pulse_sequence(5,:))
    ylabel('Gy') % in T/m

    a6 = subplot(7,1,6);
    plot(timeaxis,pulse_sequence(6,:))
    ylabel('Gz') % in T/m

    a7 = subplot(7,1,7);
    plot(timeaxis,isInKspace,'k')
    hold on
    plot(timeaxis(1,kspace_times(:)),isInKspace(1,kspace_times(:)),'ro')
    ylabel('isInKspace')
    set(gcf,'Name','Pulse Sequence')

    linkaxes([a1 a2 a3 a4 a5 a6 a7],'x');
    set(gca,'XTick',[])
end