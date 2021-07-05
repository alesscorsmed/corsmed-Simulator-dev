function [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,...
    N_TR,kspace_times,BW1,TR_times,TRmin,PulseqDuration,encodedkSpace] = ...
    SEgeneratorWithIRoption(RFmatrix1,RFmatrix2,acquiredkspace,FOV,receiver_BW,G_max,TR,...
    TE,dt,gamma,conn_localdb,experiment_id,pulseq_id,partialFourierStruct,...
    isIR,TI,IRtype,IRduration,struct_pulseq)
% This function generates Spin Echo (SE) pulse sequences. It does NOT take 
% into account GR ramps.

N_TR                = floor(TR/dt);   % Number of time-steps during the TR

rf_dur1             = RFmatrix1(1,1);
cycles1             = RFmatrix1(2,1);
angle1              = RFmatrix1(3,1);
slice_thickness1    = RFmatrix1(4,1);

Dkx                 = 1/FOV(1,1);   %k-space Dkx   ALLAXE TO
Dky                 = 1/FOV(1,2);   %k-space Dky   ALLAXE TO

% N_pulse             = kspace(1,2)*N_TR;
% G_pulse             = zeros(3,N_pulse);
% B1_pulse            = zeros(2,N_pulse);

% Acquire less k-space lines if halfFourierFactor < 1
partialFourierFactor    = partialFourierStruct.factor;
if strcmp(partialFourierStruct.type,'phaseConjugate')
    acquiredLines   = round(partialFourierFactor*acquiredkspace(1,2));
    acquiredRO      = acquiredkspace(1,1);
    
    readConjugateFactor = 1;

    if partialFourierFactor<1
        fprintf('Half-fourier is active with an acceleration factor of %.2f\n',partialFourierFactor)
        fprintf('%i lines will be acquired (out of %i)\n',acquiredLines,acquiredkspace(1,2))
    end
elseif strcmp(partialFourierStruct.type,'readConjugate')
    acquiredLines   = acquiredkspace(1,2);
    acquiredRO      = round(partialFourierFactor*acquiredkspace(1,1));
    
    readConjugateFactor = partialFourierFactor;

    if partialFourierFactor<1
        fprintf('Half-fourier is active with an acceleration factor of %.2f\n',partialFourierFactor)
        fprintf('%i samples will be acquired (out of %i)\n',acquiredRO,acquiredkspace(1,1))
    end
else
    acquiredLines   = acquiredkspace(1,2);
    acquiredRO      = acquiredkspace(1,1);
    
    readConjugateFactor = 1;
end

flag                = 0;

rf_dur2             = RFmatrix2(1,1);
cycles2             = RFmatrix2(2,1);
angle2              = RFmatrix2(3,1);
slice_thickness2    = RFmatrix2(4,1);

%% Inversion RF pulse
if isIR
    % TI counts from the end of the IR pulse till the start of the host
    % pulse sequence
    TIendPoint          = 0;

    structIR.TI         = str2num(struct_pulseq.ir_time);
    structIR.type       = struct_pulseq.ir_type;
    structIR.IRduration = str2num(struct_pulseq.ir_duration);
    structIR.angle      = 180;
    structIR.phase      = 90; % in degrees
    structIR.IRcycles   = 2; % only for sinc RF pulse

    structExper.dt      = dt;
    structExper.gamma   = gamma;

    compressIRmodule    = 1;

    % blockIR_actual_timesteps holds the "true" number of timesteps of the
    % blockIR. blockIR_timesteps holds the "effective" number of timesteps due
    % to the application of the fast algorithm

    [blockIR,blockIR_actual_timesteps,blockIR_timesteps] = ...
        pulseSequenceGenerator.generateIRblock(structIR,...
        structExper,TIendPoint,'SE',compressIRmodule,conn_localdb,...
        experiment_id,pulseq_id);
else
    blockIR_actual_timesteps    = 0;
    blockIR_timesteps           = 0;
    blockIR                     = [];
end

%% 1st Block - RF90 and Gz (slice selection)

% STEP 1 - Design RF
[RF1,BW1,rf_timesteps1]    = pulseSequenceGenerator.addRFsinc(rf_dur1,cycles1,...
    angle1,dt,gamma);

% STEP 2 - Design slice selection gradient (GZ)
Gz_magn1                 = BW1/(gamma*slice_thickness1);
if Gz_magn1>G_max
    msg = ['The strength of the slice selection gradient exceeds the ',...
            'maximum available gradient strength of the Virtual MR scanner. ',...
            'Please consider increasing the duration of the RF pulse or ',...
            'the slice thickness.'];

    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end
[GZ1,GZ_plateau1,GZ_area1] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,Gz_magn1,rf_dur1,...
    0,'GZ',dt);

% Check Error: if the length of RF pulse is equal to the length of the plateau
size_plateau1 = size(GZ_plateau1(1,2):GZ_plateau1(1,3),2);
if size_plateau1~=rf_timesteps1
    msg = ['ERROR - DO NOT IGNORE IT - The length of RF pulse is NOT',...
        ' equal to the length of the plateau'];
    disp(msg)
end

block1 = zeros(8,max(size(RF1,2),size(GZ1,2)));
block1(6,1:size(GZ1,2)) = GZ1;
block1(1,GZ_plateau1(1,2):GZ_plateau1(1,3)) = RF1(1,:);

% The phase of the first RF pulse is -90.
block1(2,GZ_plateau1(1,2):GZ_plateau1(1,3)) = 90*pi/180; % phase -90 (in rads)

disp(['Gss 90RF: ',num2str(GZ_area1)])

%% 3rd Block - RF180 and Gz (slice selection)

% STEP 1 - Design RF
[RF2,BW2,rf_timesteps2]    = pulseSequenceGenerator.addRFsinc(rf_dur2,cycles2,...
    angle2,dt,gamma);

% STEP 2 - Design slice selection gradient (GZ)
Gz_magn2                 = BW2/(gamma*slice_thickness2);
if Gz_magn2>G_max
    msg = ['The strength of the slice selection gradient exceeds the ',...
            'maximum available gradient strength of the Virtual MR scanner. ',...
            'Please consider increasing the duration of the RF pulse or ',...
            'the slice thickness.'];

    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end
[GZ2,GZ_plateau2,GZ_area2] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,Gz_magn2,rf_dur2,...
    0,'GZ',dt);

disp(['Gss 180RF: ',num2str(GZ_area2)])

% Check Error: if the length of RF pulse is equal to the length of the plateau
size_plateau2 = size(GZ_plateau2(1,2):GZ_plateau2(1,3),2);
if size_plateau2~=rf_timesteps2
    msg = ['ERROR - DO NOT IGNORE IT - The length of RF pulse is NOT',...
        ' equal to the length of the plateau'];
    disp(msg)
end

block3 = zeros(8,max(size(RF2,2),size(GZ2,2)));
block3(6,1:size(GZ2,2)) = GZ2;
block3(1,GZ_plateau2(1,2):GZ_plateau2(1,3)) = RF2(1,:);


%% 4th Block - GR_X readout

% STEP 5 - GR_X
GX_magn                 = receiver_BW/(gamma*FOV(1,1));
GX_plateau_duration     = acquiredRO/receiver_BW;  
[GX,GX_plateau,GX_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_magn,...
    GX_plateau_duration,0,'GX',dt);

block4 = zeros(8,size(GX,2));
block4(4,1:size(GX,2)) = GX;

%% 2nd Block - GR_Z ref. + GR_X pref.

% STEP 3 - GZ refocussing gradient (GZ_ref)
% We assume that the refocusing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "slice selection" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GZ_area(1,2)/2 + GZ_area(1,3))

GZ_ref_magn             = 0.95*G_max;                                                                                     
Gz_ref_plateau_tmstps   = floor((GZ_area1(1,2)/2 + GZ_area1(1,3))/...
    (GZ_ref_magn*dt));
GZ_ref_magn             = (GZ_area1(1,2)/2 + GZ_area1(1,3))/...
    (Gz_ref_plateau_tmstps*dt);  % recalculate GZ_ref_magn so as to yield 
        % equal gradient areas for the time Gz_ref_plateau_tmstps

% [GZ_ref,GZ_ref_plateau,GZ_ref_area] = ...
%     pulseSequenceGenerator.addGRtrapezoid(0,...
%     (GZ_ref_magn),Gz_ref_plateau_tmstps*dt,0,'GZ_ref',dt);
% GZ_ref = -GZ_ref;
% GZ_ref_area = -GZ_ref_area;
% 
% disp(['GssRef 90RF: ',num2str(GZ_ref_area)])

% STEP 4 - GR_X pre.
% We assume that the pre-phasing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,1)). The dephasing lobe of the frequency 
% encode (FE) gradient between the two RF pulses is now positive. This is 
% because in SE the RF does the rephasing by advancing or retarding the 
% angle of the transverse magnetization components, which continue their
% phase evolution in the same rotational sense as before.

GX_pre_magn = 0.95*G_max;
if strcmp(partialFourierStruct.type,'readConjugate')
    GxMainPart = (readConjugateFactor - 0.5) * (GX_area(1,2)/readConjugateFactor);
else
    GxMainPart = GX_area(1,2)/2;
end
GX_pre_plateau_tmstps = ceil((GxMainPart + GX_area(1,1))/(GX_pre_magn*dt));
GX_pre_magn = (GxMainPart + GX_area(1,1))/(GX_pre_plateau_tmstps*dt);

[GX_pre,GX_pre_plateau,GX_pre_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_pre_magn,...
    GX_pre_plateau_tmstps*dt,0,'GX_pre',dt);

%% TEMP

GZ2_ref_magn             = 0.95*G_max;                                                                                    
Gz2_ref_plateau_tmstps   = floor((GZ_area2(1,2)/2 + GZ_area2(1,3))/...
    (GZ2_ref_magn*dt));
G2Z_ref_magn             = (GZ_area2(1,2)/2 + GZ_area2(1,3))/...
    (Gz2_ref_plateau_tmstps*dt);  % recalculate GZ_ref_magn so as to yield 
        % equal gradient areas for the time Gz_ref_plateau_tmstps

[GZ2_ref,GZ2_ref_plateau,GZ2_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    GZ2_ref_magn,Gz2_ref_plateau_tmstps*dt,0,'GZ2_ref',dt);

disp(['GssRef 180RF: ',num2str(GZ2_ref_area)])

[GZ_ref,GZ_ref_plateau,GZ_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    GZ2_ref_magn,Gz2_ref_plateau_tmstps*dt,0,'GZ_ref',dt);
GZ_ref = -GZ_ref;
GZ_ref_area = -GZ_ref_area;

disp(['GssRef 90RF: ',num2str(GZ_ref_area)])

%% Calculate the min TE

block2_minDuration = (max(size(GZ_ref,2)+size(GZ2_ref,2),size(GX_pre,2)))*dt;
TEmin = 2*((size(block1,2)/2 + size(block3,2)/2)*dt + block2_minDuration);

% If there is an overlap between the second half of the second RF and the
% first (till TE) part of readout
if TEmin/2 < (size(block3,2)/2 + (1 - 0.5/readConjugateFactor)*size(block4,2))*dt
    TEminProposed = 2*(size(block3,2)/2 + (1 - 0.5/readConjugateFactor)*size(block4,2))*dt;
    msg = ['The selected TE is smaller than the minimum available TE (',...
        num2str(TEminProposed),'sec) for the current configuration. Consider',...
        ' increasing the TE or decreasing the duration of the readout.'];
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end

% Max. area of GY_PE
GY_ramp         = 0;
GY_max          = 0.95*G_max;
ky(1,1)         = ((acquiredkspace(1,2)-1)/2)*Dky;
GY_area(1,1) 	= ky(1,1)/gamma;

% Crusher gradients will be applied right before and after the 180 RF
% pulse. The crushers will be equal to the GZ2_ref
GYcrusher = GZ2_ref;
block2_minDuration = block2_minDuration - size(GYcrusher,2)*dt;

while GY_area(1,1)/block2_minDuration > GY_max
    
    block2_minDuration  = block2_minDuration + dt;
    TEmin               = TEmin + dt;
    flag                = 1;

end

if TE < TEmin
    
    if flag == 1
        
        msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or the maximum gradient strength (if applicable).'];        
        
    else
    
%         msg = ['The selected TE is smaller than the minimum available TE (',...
%             num2str(TEmin),'sec) for the current configuration. Consider',...
%             ' decreasing the duration of the RF pulse or the duration of',...
%             ' the readout.'];
        msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or decreasing the duration of the RF pulse or the duration of',...
            ' the readout (this option may not be available).'];
        
    end
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end

%% 2nd Block - PE
% TE is equal to 2*(block1/2+block2+block3/2) => 
% block2 = (TE - block1 - block3)/2
block2_Duration = (TE - (size(block1,2) + size(block3,2))*dt)/2;

% In this block, crusher gradient will be applied right before the 180 RF
% pulse. The crusher will be equal to the GZ2_ref
block2_GY_Duration = block2_Duration - size(GYcrusher,2)*dt;

if block2_GY_Duration < block2_minDuration
    
    TEmin = 2*(block2_minDuration + size(GYcrusher,2)*dt) + ...
        (size(block1,2) + size(block3,2))*dt;
    
    TEmin = ceil(TEmin*1000)/1000;
    
    msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or the maximum gradient strength (if applicable).'];  
        
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
        
end

% STEP 6 - PE
if GY_area(1,1)/block2_GY_Duration <= GY_max
    
    if GY_area(1,1)/GY_max < dt
        
        msg = ['The selected dt is larger than the needed dt (',...
            num2str(GY_area(1,1)/GY_max),' sec) for running the max. PE ',...
            ' gradient. Consider decreasing the dt.'];
        
        messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
        
    end

    GY_magn = GY_area(1,1)/block2_GY_Duration;
    
    [GY,GY_plateau,GY_area_1] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,block2_GY_Duration,0,'GY1',dt);
    

    %------CONTROL-------%
    % Check if the area of the Gy_max is finally equal to the wanted area
    % of Gy_area(1,n). I am using the str2num(num2str()) convertion since
    % isequal() doesn't work
    GY_area_total = sum(GY_area_1,2);
    if ~isequal(str2num(num2str(GY_area(1,1))),...
            str2num(num2str(GY_area_total))) %#ok
        disp('Gy: PROBLEM WITH THE AREA OF maxGy');
    end
    %------END CONTROL-------%    
    
end

block2 = zeros(8,floor(block2_Duration/dt));
block2(6,1:size(GZ_ref,2)) = GZ_ref;
block2(4,end-size(GX_pre,2)+1:end) = GX_pre;
block2(5,1:size(GY,2)) = GY;
block2(5,end-size(GYcrusher,2)+1:end) = GYcrusher;
block2(6,end-size(GZ2_ref,2)+1:end) = GZ2_ref;

%% BUILD TR
TR_actual   = blockIR_timesteps + size(block1,2)/2+round(TE/dt)+...
    round(size(block4,2)*(0.5/readConjugateFactor));
TRmin       = blockIR_actual_timesteps + size(block1,2)/2+round(TE/dt)+...
    round(size(block4,2)*(0.5/readConjugateFactor));

if N_TR<TRmin
    
    msg = ['The selected TR is lower that the minimum available TR (',...
            num2str(TRmin*dt),' sec) for the current configuration. Consider ',...
            'increasing the TR.'];
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end

pulse_sequence_TR = zeros(8,TR_actual);
current_timestep = 0;
if isIR
    pulse_sequence_TR(:,current_timestep+1:current_timestep+size(blockIR,2)) = blockIR;
    current_timestep = current_timestep + size(blockIR,2);
end
pulse_sequence_TR(:,current_timestep+1:current_timestep+size(block1,2)) = block1;
current_timestep = current_timestep + size(block1,2);
pulse_sequence_TR(:,current_timestep+1:...
    current_timestep+size(block2,2)) = block2;
current_timestep = current_timestep + size(block2,2);
pulse_sequence_TR(:,current_timestep+1:current_timestep+size(block3,2)) = block3;
current_timestep = current_timestep+size(block3,2);
pulse_sequence_TR(:,TR_actual-size(block4,2)+1:TR_actual) = block4;

% pulse_sequence_TR(:,1:size(block1,2)) = block1;
% 
% pulse_sequence_TR(:,size(block1,2)+1:...
%     size(block1,2)+size(block2,2)) = block2;
% 
% pulse_sequence_TR(:,size(block1,2)+...
%     size(block2,2)+1:size(block1,2)+size(block2,2)+size(block3,2)) = block3;
% pulse_sequence_TR(:,TRmin-size(block4,2)+1:TRmin) = block4;

% pulse_sequence_TR(7,:) = 1:size(pulse_sequence_TR,2);
% pulse_sequence_TR(8,:) = 1:size(pulse_sequence_TR,2);

pulse_sequence_TR(7,size(blockIR,2)+1:end) = size(blockIR,2)+1:size(pulse_sequence_TR,2);
pulse_sequence_TR(8,size(blockIR,2)+1:end) = size(blockIR,2)+1:size(pulse_sequence_TR,2);

%%
pulse_sequence_TR(6,current_timestep+1:...
    current_timestep+size(GZ2_ref,2)) = GZ2_ref;

% Add the crusher gradient right after the 180 RF pulse
pulse_sequence_TR(5,current_timestep+1:...
    current_timestep+size(GYcrusher,2)) = GYcrusher;

%%

% Add two more extra points. The first holds the accumulated effect of the
% pulse sequence till the end of the TR (minus one dt) and the second holds
% the temporal effect of the soft crushers at the end of the TR
pulse_sequence_TR           = [pulse_sequence_TR,[zeros(6,1);TR_actual+1;N_TR-1]];
pulse_sequence_TR           = [pulse_sequence_TR,[zeros(6,1);N_TR;N_TR]];

soft_crushers_TR            = zeros(1,size(pulse_sequence_TR,2));
soft_crushers_TR(1,size(pulse_sequence_TR,2))   = 1;
if isIR
    % crusher at the end of the IR block is needed
    soft_crushers_TR(1,size(blockIR,2)) = 1;
end

pulse_sequence              = repmat(pulse_sequence_TR,1,acquiredLines);
soft_crushers               = repmat(soft_crushers_TR,1,acquiredLines);

isInKspace                  = zeros(1,size(pulse_sequence,2));

% % Delete all the other PE steps. You will create them later in this script
% pulse_sequence(5,size(pulse_sequence_TR,2)+1:end) = ...
%     zeros(1,size(size(pulse_sequence_TR,2)+1:size(pulse_sequence,2),2));

kspace_times                = zeros(acquiredLines,2);
TR_times                    = zeros(acquiredLines,2);

temp1 = 0;
for i=1:acquiredLines
    
    % Due to the implementation of the fast algorithm at the end of each
    % TR, the length of each TR is not N_TR but TRmin+2
    kspace_times(i,1) = ((i-1)*(TR_actual+2)) + TR_actual - size(block4,2) + 1;
    kspace_times(i,2) = kspace_times(i,1) + (GX_plateau(1,3)-GX_plateau(1,2));
    
    TR_times(i,1) = (i-1)*(TR_actual+2)+1;
    TR_times(i,2) = i*(TR_actual+2);
    
    isInKspace(1,kspace_times(i,1):kspace_times(i,2)) = temp1 + [1:GX_plateau(1,3)];
    temp1 = i*GX_plateau(1,3);
end

for n=2:acquiredLines
    
    pulse_sequence(7,(n-1)*size(pulse_sequence_TR,2)+1:...
        n*size(pulse_sequence_TR,2)) = (n-1)*N_TR+pulse_sequence_TR(7,:);
    pulse_sequence(8,(n-1)*size(pulse_sequence_TR,2)+1:...
        n*size(pulse_sequence_TR,2)) = (n-1)*N_TR+pulse_sequence_TR(8,:);
    
    if ((acquiredkspace(1,2)-1)/2-(n-1))==0   % CONTROL: In case we have an odd 
        % number of kspace lines, the Gy for the line that passes through 
        % zero must have the following characteristics
        ky(1,n) = 0;
        GY_area(1,n) = 0; 
        GY_magn = 0;
    else
        ky(1,n)=((acquiredkspace(1,2)-1)/2-(n-1))*Dky;
        GY_area(1,n) = -ky(1,n)/gamma;  
        GY_magn = GY_area(1,n)/block2_GY_Duration;
    end
    
    [GY,GY_plateau,GY_area_2] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,block2_GY_Duration,0,['GY',num2str(n)],dt);
                                           
    %------CONTROL-------%
    % Check if the area of the Gy_max is finally equal to the wanted area
    % of Gy_area(1,n). I am using the str2num(num2str()) convertion since
    % isequal() doesn't work
    GY_area_total = sum(GY_area_2,2);
    if ~isequal(str2num(num2str(GY_area(1,n))),...
            str2num(num2str(GY_area_total))) %#ok
        disp(['Gy: PROBLEM WITH THE AREA OF GY',num2str(n)]);
    end
    %------END CONTROL-------%
    
    pulse_sequence(5,TR_times(n-1,2)+size(blockIR,2)+size(block1,2)+1:...
        TR_times(n-1,2)+size(blockIR,2)+size(block1,2)+size(GY,2)) = GY;
    
end

N_pulse = size(pulse_sequence,2);

%% Calculate extra parameters
encodedkSpace.lines         = acquiredLines;
encodedkSpace.columns       = round(GX_plateau_duration/dt);

PulseqDuration              = (TR + blockIR_actual_timesteps*dt)*acquiredLines;
disp(['The total duration of the pulse sequence is ',num2str(PulseqDuration),'sec'])

%% FINAL CHECK
if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
        max(pulse_sequence(6,:))>=G_max
    exception = MException('GrMagn:GreaterThanMax',...
        'The gradients magnitude exceeds the maximum gradient strength');
    throw(exception)
end