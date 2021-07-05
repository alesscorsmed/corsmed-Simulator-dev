%% INITIAL CONFIGURATION
% Model-specific parameters
tic
model_desired_sizeofelement = [0.0005,0.0005,0.0005];
model_desired_dim = [0.3, 0.2, 0.3];

gamma           = 42.56*10^6;   % Hz/T
dt              = 10^-6;        % sec
sliceThickness  = 0.006;        % desired slice thickness (in m)
f0              = 64000000;     % in Hz

fast            = 0;
save            = 0;
plotON          = 1;

RFmatrix1        = zeros(6,1);
RFmatrix1(1,:)   = 0.003;            % RF duration 
RFmatrix1(2,:)   = 2;                % RF sinc cycles
RFmatrix1(3,:)   = 90;               % RF rotation angle
RFmatrix1(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
RFmatrix1(5,:)   = [1];              % starting timestep of RF
RFmatrix1(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis

RFmatrix2        = zeros(6,1);
RFmatrix2(1,:)   = 0.003;            % RF duration 
RFmatrix2(2,:)   = 2;                % RF sinc cycles
RFmatrix2(3,:)   = 180;              % RF rotation angle
RFmatrix2(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
RFmatrix2(5,:)   = [1];              % starting timestep of RF
RFmatrix2(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis

receiverBW      = 200*10^3;  % in Hz (this is a feature of the hardware)
G_max           = 40*10^-3; % in T/m

kspace(1,1) = 128; %defines the size of k-space along x direction  %128x128
kspace(1,2) = 128; %defines the size of k-space along y direction

FOV(1,1)    = 0.4;
FOV(1,2)    = 0.4;

TR          = 4;     % sec
TE          = 0.06;

acquiredkspace      = kspace;
acquiredFOV         = FOV;

partialFourierStruct.type   = 'readConjugate';
partialFourierStruct.factor = 0.8;

receiver_BW_star = 1/dt; % this is the BW which depends on the time step 
                        % which is chosen for the simulation. It is needed
                        % to calculate the factor_BW by which the kspace(1,1)
                        % will be multiplied. This product will define the
                        % points that will be collected from the signal
                        % during the acquisition.

[pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
    kspace_times,BW_rf,TR_times,TRmin,PulseqDuration,encodedkSpace] = ...
    pulseSequenceGenerator.SEgenerator(RFmatrix1,RFmatrix2,kspace,FOV,...
    receiverBW,G_max,TR,TE,dt,gamma,[],[],[],partialFourierStruct);
  
%% ALARM - CHECK IF SPURIOUS ECHOES WILL APPEAR
% PHASE DIFFERENCE OF CONTIGUOUS ELEMENTS ALARM
elementsize=model_desired_sizeofelement;
alarm=[1,1,1];
alarm_count = 0;
while alarm(1,1)==1||alarm(1,2)==1||alarm(1,3)==1    
    [elementsize,alarm] = ...
        pulseSequenceGenerator.testAccumPhaseGRE(N_pulse,...
        N_TR,dt,pulse_sequence(4:6,:),elementsize);
    alarm_count = alarm_count + 1;
end

disp(['The alarm was enabled ',num2str(alarm_count - 1),' times.'])
disp(['The elementsize should be ',num2str(elementsize(1,1)),' x ',...
    num2str(elementsize(1,2)),' x ',num2str(elementsize(1,3)),' meters (m)'])

%% FIX THE PULSE SEQUENCE TO WORK PROPERLY WITH coreMRI
info = pulseSequenceGenerator.createInfoStruct_SE(kspace,FOV,dt,...
    receiverBW,TR,TE,kspace_times,f0,sliceThickness,acquiredkspace,...
    acquiredFOV,encodedkSpace);

%% ADD FAST ALGORITHM
if fast == 1
    [pulse_sequence_new,kspace_times_new,isInKspace_new,soft_crushers_new] = ...
        pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
        info.pulseSequence.kspace,soft_crushers);
    
    info_new = pulseSequenceGenerator.createInfoStruct_SE(kspace,FOV,dt,...
        receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,acquiredkspace,...
        acquiredFOV,encodedkSpace);
end

%% SAVE PULSE SEQUENCE
if save == 1
    
    pulseSequenceGenerator.savePulseSequence(pulse_sequence,...
        soft_crushers,isInKspace,N_pulse,info,gamma,dt,'test')
    
    if fast == 1
        pulseSequenceGenerator.savePulseSequence(pulse_sequence_new,...
            soft_crushers_new,isInKspace_new,size(pulse_sequence_new,2),...
            info_new,gamma,dt,'test_fast')
    end
    
end

%% PLOT
if plotON == 1
    timeaxis = [1:N_pulse]*dt;

    a1 = subplot(7,1,1);
    plot(timeaxis,pulse_sequence(1,:),'r')
    ylabel('RF magnitude') % in Tesla

    a2 = subplot(7,1,2);
    plot(timeaxis,pulse_sequence(2,:),'r')
    ylabel('RF phase') % in radians

    a3 = subplot(7,1,3);
    plot(timeaxis,pulse_sequence(3,:),'r')
    ylabel('RF Freq.') % in Hz

    a4 = subplot(7,1,4);
    plot(timeaxis,pulse_sequence(4,:))
    ylabel('Gx') % in T/m

    a5 = subplot(7,1,5);
    plot(timeaxis,pulse_sequence(5,:))
    ylabel('Gy') % in T/m

    a6 = subplot(7,1,6);
    plot(timeaxis,pulse_sequence(6,:))
    ylabel('Gz') % in T/m

    a7 = subplot(7,1,7);
    plot(timeaxis,isInKspace,'k')
    ylabel('isInKspace')
    set(gcf,'Name','Pulse Sequence')

    linkaxes([a1 a2 a3 a4 a5 a6 a7],'x');
    set(gca,'XTick',[])
end
toc