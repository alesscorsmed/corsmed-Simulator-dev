function [pulse_sequence,isInKspace,soft_crushers,N_pulse,...
    kspace_times,BW2,PulseqDuration,ESP,encodedkSpace] = ...
    EPI_SEgenerator(RFmatrix1,RFmatrix2,acquiredkspace,FOV,receiver_BW,...
    G_max,dt,TE,gamma,conn_localdb,experiment_id,pulseq_id,partialFourierStruct)
% This function generates SE-EPI pulse sequences. It does NOT take into
% account GR ramps. ESP stands for Echo Spacing. ESP defines the duration
% of block3.
% In this file, the part of the pulse sequence that starts at
% the onset of the 180 RF pulse is first structured and then the first part
% of the pulse sequence is added at the front of the pulse sequence design.
% This pulse sequence consists of the following 5 blocks:
% block_pre1: it holds the RF90 and Gz (slice selection)
% block_pre2: it holds the GR_Z ref. (for the 90 RF pulse), the GR_X pref.,
% the GR_Y pref and the GR_Z pref. (for the 180 RF pulse)
% block_1: it holds the RF (180) and Gz (for the 180 RF pulse)
% block_2: it holds the GR_Z ref. (for the 180 RF pulse)
% block_3: it holds the GR_X readout and the GR_Y blip. block 3 is repeated
% for every kspace line. 

rf_dur1          = RFmatrix1(1,1);
cycles1          = RFmatrix1(2,1);
angle1           = RFmatrix1(3,1);
slice_thickness1 = RFmatrix1(4,1);

rf_dur2          = RFmatrix2(1,1);
cycles2          = RFmatrix2(2,1);
angle2           = RFmatrix2(3,1);
slice_thickness2 = RFmatrix2(4,1);

Dkx             = 1/FOV(1,1);   %k-space Dkx   ALLAXE TO
Dky             = 1/FOV(1,2);   %k-space Dky   ALLAXE TO

flag            = 0;

% Acquire less k-space lines if halfFourierFactor < 1
partialFourierFactor    = partialFourierStruct.factor;
acquiredLines           = round(partialFourierFactor*acquiredkspace(1,2));

if partialFourierFactor<1
    fprintf('Half-fourier is active with an acceleration factor of %.2f\n',partialFourierFactor)
    fprintf('%i line will be acquired (out of %i)\n',acquiredLines,acquiredkspace(1,2))
end

GY_blip_area    = Dky/gamma;
GY_App          = ((acquiredkspace(1,2)-1)/2)*GY_blip_area;

%% 1st Block - RF (180) and Gz (slice selection)

% STEP 1 - Design RF
[RF2,BW2,rf2_timesteps]    = pulseSequenceGenerator.addRFsinc(rf_dur2,cycles2,...
    angle2,dt,gamma);

% STEP 2 - Design slice selection gradient (GZ)
Gz_magn                 = BW2/(gamma*slice_thickness2);
if Gz_magn>G_max
    msg = ['The strength of the slice selection gradient exceeds the ',...
            'maximum available gradient strength of the Virtual MR scanner. ',...
            'Please consider increasing the duration of the RF pulse or ',...
            'the slice thickness.'];

    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end
[GZ2,GZ2_plateau,GZ2_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,Gz_magn,rf_dur2,...
    0,'GZ',dt);

disp(['Gss 180RF: ',num2str(GZ2_area)])

% Check Error: if the length of RF pulse is equal to the length of the plateau
size_plateau = size(GZ2_plateau(1,2):GZ2_plateau(1,3),2);
if size_plateau~=rf2_timesteps
    msg = ['ERROR - DO NOT IGNORE IT - The length of RF pulse is NOT',...
        ' equal to the length of the plateau'];
    disp(msg)
end

block1 = zeros(8,max(size(RF2,2),size(GZ2,2)));
block1(6,1:size(GZ2,2)) = GZ2;
block1(1,GZ2_plateau(1,2):GZ2_plateau(1,3)) = RF2(1,:);


%% 3rd Block - GR_X readout + GR_Y blip

% STEP 6 - GR_X
GX_magn                 = receiver_BW/(gamma*FOV(1,1));
GX_plateau_duration     = acquiredkspace(1,1)/receiver_BW;  
[GX,GX_plateau,GX_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_magn,...
    GX_plateau_duration,0,'GX',dt);

% STEP 7 - GR_Y blip
GYblip_magn                 = GY_blip_area/dt;
GYblip_plateau_duration     = dt;  
[GYblip,GYblip_plateau,GYblip_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GYblip_magn,...
    GYblip_plateau_duration,0,'GYblip',dt);

block3 = zeros(8,size(GX,2)+1);
block3(4,1:size(GX,2)) = GX;
block3(5,end) = GYblip;

ESP = size(block3,2)*dt;

%% 2nd Block - GR_Z ref.

% GZ refocussing gradient (GZ_ref) (180 RF pulse)
% We assume that the refocusing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "slice selection" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GZ_area(1,2)/2 + GZ_area(1,3))

GZ2_ref_magn             = 0.95*G_max;                                                                                    
GZ2_ref_plateau_tmstps   = floor((GZ2_area(1,2)/2 + GZ2_area(1,3))/...
    (GZ2_ref_magn*dt));
GZ2_ref_magn             = (GZ2_area(1,2)/2 + GZ2_area(1,3))/...
    (GZ2_ref_plateau_tmstps*dt);  % recalculate GZ_ref_magn so as to yield 
        % equal gradient areas for the time Gz_ref_plateau_tmstps

[GZ2_ref,GZ2_ref_plateau,GZ2_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    GZ2_ref_magn,GZ2_ref_plateau_tmstps*dt,0,'GZ2_ref',dt);

disp(['GssRef 180RF: ',num2str(GZ2_ref_area)])

% Calculate the min duration of block2
block2_minDuration = size(GZ2_ref,2)*dt;

block2 = zeros(8,floor(block2_minDuration/dt));
block2(6,1:size(GZ2_ref,2)) = GZ2_ref;

% Check the minimum echo time
% The minimum echo time (TE) is defined as two-times the time between the
% center of the 180 RF pulse (this is why we devide the size of block1
% by two) and the acquisition of the center of the kspace.
TEmin = 2*...
    (round(size(block1,2)/2)+size(block2,2)+(acquiredkspace(1,2)/2)*size(block3,2))*dt;

if TE < TEmin
        
    msg = ['The selected TE is smaller than the minimum available TE (',...
        num2str(TEmin),'sec) for the current configuration. Consider',...
        ' increasing the TE or decreasing the number of kspace lines.'];
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',...
        experiment_id,pulseq_id);
    
end

%% Block pre 1
% STEP 1 - Design RF
[RF1,BW1,rf_timesteps1]    = pulseSequenceGenerator.addRFsinc(rf_dur1,cycles1,...
    angle1,dt,gamma);

% STEP 2 - Design slice selection gradient (GZ)
Gz1_magn1                 = BW1/(gamma*slice_thickness1);
[GZ1,GZ1_plateau,GZ1_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,Gz1_magn1,rf_dur1,...
    0,'GZ',dt);

% Check Error: if the length of RF pulse is equal to the length of the plateau
size_plateau1 = size(GZ1_plateau(1,2):GZ1_plateau(1,3),2);
if size_plateau1~=rf_timesteps1
    msg = ['ERROR - DO NOT IGNORE IT - The length of RF pulse is NOT',...
        ' equal to the length of the plateau'];
    disp(msg)
end

block_pre1 = zeros(8,max(size(RF1,2),size(GZ1,2)));
block_pre1(6,1:size(GZ1,2)) = GZ1;
block_pre1(1,GZ1_plateau(1,2):GZ1_plateau(1,3)) = RF1(1,:);

% The phase of the first RF pulse is -90.
block_pre1(2,GZ1_plateau(1,2):GZ1_plateau(1,3)) = 90*pi/180; % phase -90 (in rads)

disp(['Gss 90RF: ',num2str(GZ1_area)])

%% Block pre 2:
% GZ refocussing gradient (GZ_ref) for the 90 RF pulse
% We assume that the refocusing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "slice selection" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GZ_area(1,2)/2 + GZ_area(1,3))

GZ1_ref_magn             = 0.95*G_max;                                                                                     
Gz1_ref_plateau_tmstps   = floor((GZ1_area(1,2)/2 + GZ1_area(1,3))/...
    (GZ1_ref_magn*dt));
GZ1_ref_magn             = (GZ1_area(1,2)/2 + GZ1_area(1,3))/...
    (Gz1_ref_plateau_tmstps*dt);  % recalculate GZ_ref_magn so as to yield 
        % equal gradient areas for the time Gz_ref_plateau_tmstps
        
[GZ1_ref,GZ1_ref_plateau,GZ1_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    -GZ1_ref_magn,Gz1_ref_plateau_tmstps*dt,0,'GZ1_ref',dt);

disp(['GssRef 90RF: ',num2str(GZ2_ref_area)])

% GR_X pre.
% We assume that the pre-phasing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,1))

GX_pre_magn = 0.95*G_max;                                                 
GX_pre_plateau_tmstps = floor((GX_area(1,2)/2 + GX_area(1,1))/(GX_pre_magn*dt));
GX_pre_magn = (GX_area(1,2)/2 + GX_area(1,1))/(GX_pre_plateau_tmstps*dt);

[GX_pre,GX_pre_plateau,GX_pre_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_pre_magn,...
    GX_pre_plateau_tmstps*dt,0,'GX_pre',dt);

% GR_Y pre.
% We assume that the pre-phasing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,1))
% GY_App
GY_pre_magn = 0.95*G_max;                                                 
GY_pre_plateau_tmstps = floor(abs(GY_App)/(GY_pre_magn*dt));
GY_pre_magn = (GY_App)/(GY_pre_plateau_tmstps*dt);

[GY_pre,GY_pre_plateau,GY_pre_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GY_pre_magn,...
    GY_pre_plateau_tmstps*dt,0,'GY_pre',dt);

block_pre2_minDuration = max([size(GY_pre,2),size(GX_pre,2),...
    (size(GZ1_ref,2)+size(GZ2_ref,2))])*dt;

part1duration = block_pre2_minDuration + round(size(block_pre1,2)/2)*dt + ...
    round(size(block1,2)/2)*dt;

if part1duration < TE/2
    
    part1duration = TE/2;
    block_pre2_Duration = part1duration - round(size(block_pre1,2)/2)*dt - ...
        round(size(block1,2)/2)*dt;
    
end

block_pre2 = zeros(8,floor(block_pre2_Duration/dt));
block_pre2(6,1:size(GZ1_ref,2)) = GZ1_ref;
block_pre2(6,end-size(GZ2_ref,2)+1:end) = GZ2_ref;
block_pre2(4,end-size(GX_pre,2)+1:end) = GX_pre;
block_pre2(5,end-size(GY_pre,2)+1:end) = GY_pre;

%% BUILD THE PULSE SEQUENCE
N_pulse = size(block_pre1,2)+size(block_pre2,2)+size(block1,2)+...
    size(block2,2)+acquiredLines*size(block3,2);

pulse_sequence = zeros(8,N_pulse);
pulse_sequence(:,1:size(block_pre1,2)) = block_pre1;
pulse_sequence(:,size(block_pre1,2)+1:...
    size(block_pre1,2)+size(block_pre2,2)) = block_pre2;
pulse_sequence(:,size(block_pre1,2)+size(block_pre2,2)+1:...
    size(block_pre1,2)+size(block_pre2,2)+size(block1,2)) = block1;
pulse_sequence(:,size(block_pre1,2)+size(block_pre2,2)+size(block1,2)+1:...
    size(block_pre1,2)+size(block_pre2,2)+size(block1,2)+size(block2,2)) = block2;

isInKspace = zeros(1,N_pulse);

temp1 = 0;
for i=1:acquiredLines
    startTime   = size(block_pre1,2)+size(block_pre2,2)+size(block1,2)+size(block2,2)+(i-1)*size(block3,2)+1;
    endTime     = size(block_pre1,2)+size(block_pre2,2)+size(block1,2)+size(block2,2)+i*size(block3,2);
    
    pulse_sequence(:,startTime:endTime) = block3;
    pulse_sequence(4,startTime:endTime) = (-1).^(i-1)*pulse_sequence(4,startTime:endTime);
    kspace_times(i,1) = startTime;
    kspace_times(i,2) = endTime - 1;
    
    isInKspace(1,kspace_times(i,1):kspace_times(i,2)) = temp1 + [1:size([kspace_times(i,1):kspace_times(i,2)],2)];
    temp1 = i*size([kspace_times(i,1):kspace_times(i,2)],2);
end

soft_crushers = zeros(1,N_pulse);

%% Calculate extra parameters
encodedkSpace.lines         = acquiredLines;
encodedkSpace.columns       = round(GX_plateau_duration/dt);
encodedkSpace.contrasts     = 1;  % since we acquire only one echoe per TR

PulseqDuration = N_pulse*dt;
disp(['The total duration of the pulse sequence is ',num2str(PulseqDuration),'sec'])

% %% FINAL CHECK
% if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
%         max(pulse_sequence(6,:))>=G_max
%     exception = MException('GrMagn:GreaterThanMax',...
%         'The gradients magnitude exceeds the maximum gradient strength');
%     throw(exception)
% end