function [pulse_sequence,isInKspace,soft_crushers,N_pulse,...
    kspace_times,BW,PulseqDuration,ESP,encodedkSpace] = ...
    EPIgenerator(RFmatrix,acquiredkspace,FOV,receiver_BW,G_max,dt,gamma,...
    conn_localdb,experiment_id,pulseq_id,partialFourierStruct)
% This function generates EPI pulse sequences. It does NOT take into
% account GR ramps. ESP stands for Echo Spacing. ESP defines the duration
% of block3.

rf_dur          = RFmatrix(1,1);
cycles          = RFmatrix(2,1);
angle           = RFmatrix(3,1);
slice_thickness = RFmatrix(4,1);

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
GY_App          = -((acquiredkspace(1,2)-1)/2)*GY_blip_area;

%% 1st Block - RF and Gz (slice selection)

% STEP 1 - Design RF
[RF,BW,rf_timesteps]    = pulseSequenceGenerator.addRFsinc(rf_dur,cycles,...
    angle,dt,gamma);

% STEP 2 - Design slice selection gradient (GZ)
Gz_magn                 = BW/(gamma*slice_thickness);
if Gz_magn>G_max
    msg = ['The strength of the slice selection gradient exceeds the ',...
            'maximum available gradient strength of the Virtual MR scanner. ',...
            'Please consider increasing the duration of the RF pulse or ',...
            'the slice thickness.'];

    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end
[GZ,GZ_plateau,GZ_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,Gz_magn,rf_dur,...
    0,'GZ',dt);

% Check Error: if the length of RF pulse is equal to the length of the plateau
size_plateau = size(GZ_plateau(1,2):GZ_plateau(1,3),2);
if size_plateau~=rf_timesteps
    msg = ['ERROR - DO NOT IGNORE IT - The length of RF pulse is NOT',...
        ' equal to the length of the plateau'];
    disp(msg)
end

block1 = zeros(8,max(size(RF,2),size(GZ,2)));
block1(6,1:size(GZ,2)) = GZ;
block1(1,GZ_plateau(1,2):GZ_plateau(1,3)) = RF(1,:);


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

%% 2nd Block - GR_Z ref. + GR_X pref. + GR_Y pref.

% STEP 3 - GZ refocussing gradient (GZ_ref)
% We assume that the refocusing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "slice selection" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GZ_area(1,2)/2 + GZ_area(1,3))

GZ_ref_magn             = 0.95*G_max;                                                                                        
Gz_ref_plateau_tmstps   = floor((GZ_area(1,2)/2 + GZ_area(1,3))/...
    (GZ_ref_magn*dt));
GZ_ref_magn             = (GZ_area(1,2)/2 + GZ_area(1,3))/...
    (Gz_ref_plateau_tmstps*dt);  % recalculate GZ_ref_magn so as to yield 
        % equal gradient areas for the time Gz_ref_plateau_tmstps

[GZ_ref,GZ_ref_plateau,GZ_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    (-1*GZ_ref_magn),Gz_ref_plateau_tmstps*dt,0,'GZ_ref',dt);

% STEP 4 - GR_X pre.
% We assume that the pre-phasing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,1))

GX_pre_magn = 0.95*G_max;                                                 
GX_pre_plateau_tmstps = floor((GX_area(1,2)/2 + GX_area(1,1))/(GX_pre_magn*dt));
GX_pre_magn = (GX_area(1,2)/2 + GX_area(1,1))/(GX_pre_plateau_tmstps*dt);
GX_pre_magn = GX_pre_magn*(-1);

[GX_pre,GX_pre_plateau,GX_pre_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_pre_magn,...
    GX_pre_plateau_tmstps*dt,0,'GX_pre',dt);

% STEP 5 - GR_Y pre.
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

% Calculate the min duration of block2
block2_minDuration = (max([size(GZ_ref,2),size(GX_pre,2),size(GY_pre,2)]))*dt;

block2 = zeros(8,floor(block2_minDuration/dt));
block2(6,1:size(GZ_ref,2)) = GZ_ref;
block2(4,end-size(GX_pre,2)+1:end) = GX_pre;
block2(5,end-size(GY_pre,2)+1:end) = GY_pre;

%% BUILD THE PULSE SEQUENCE
N_pulse = size(block1,2)+size(block2,2)+acquiredLines*size(block3,2);

pulse_sequence = zeros(8,N_pulse);
pulse_sequence(:,1:size(block1,2)) = block1;
pulse_sequence(:,size(block1,2)+1:...
    size(block1,2)+size(block2,2)) = block2;

isInKspace = zeros(1,N_pulse);

temp1 = 0;
for i=1:acquiredLines
    startTime   = size(block1,2)+size(block2,2)+(i-1)*size(block3,2)+1;
    endTime     = size(block1,2)+size(block2,2)+i*size(block3,2);
    
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