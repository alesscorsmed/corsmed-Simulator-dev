function [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,...
    N_TR,kspace_times,BW,TR_times,TRmin,PulseqDuration,encodedkSpace] = ...
    MPRAGEgenerator_GRE3D(RFmatrix,acquiredkspace,FOV,receiver_BW,G_max,TR,...
    TE,dt,gamma,MPRAGETR,struct_pulseq,segments,conn_localdb,...
    experiment_id,pulseq_id)

N_TR            = floor(TR/dt);
N_MPRAGETR      = floor(MPRAGETR/dt);

G_TR            = zeros(3,N_TR);
B1_TR           = zeros(2,N_TR);

rf_dur          = RFmatrix(1,1);
cycles          = RFmatrix(2,1);
angle           = RFmatrix(3,1);
slabThickness   = RFmatrix(4,1);

Dkx             = 1/FOV(1,1);           % k-space Dkx
Dky             = 1/FOV(1,2);           % k-space Dky
Dkz             = 1/slabThickness;      % k-space Dkz

acquiredLines   = acquiredkspace(1,2);
acquiredRO      = acquiredkspace(1,1);

% N_pulse         = acquiredkspace(1,2)*N_TR*acquiredkspace(1,3);
% G_pulse         = zeros(3,N_pulse);
% B1_pulse        = zeros(2,N_pulse);

flag            = 0;

%% 1st Block - RF and Gz (slice selection)

% STEP 1 - Design RF
[RF,BW,rf_timesteps]    = pulseSequenceGenerator.addRFsinc(rf_dur,cycles,...
    angle,dt,gamma);

% STEP 2 - Design slice selection gradient (GZ)
Gz_magn                 = BW/(gamma*slabThickness);
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


%% 3rd Block - GR_X readout

% STEP 5 - GR_X
GX_magn                 = receiver_BW/(gamma*FOV(1,1));
GX_plateau_duration     = acquiredkspace(1,1)/receiver_BW;  
[GX,GX_plateau,GX_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_magn,...
    GX_plateau_duration,0,'GX',dt);

block3 = zeros(8,size(GX,2));
block3(4,1:size(GX,2)) = GX;

%% 2nd Block - GR_Z ref. + GR_X pref. + slice encoding gradient

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

%% Calculate the min TE

% Check fig.11.30 at pg. 435 of Handbook of MRI pulse sequences
minimizeTE = 0;

block2_minDuration = (max(size(GZ_ref,2),size(GX_pre,2)))*dt;
TEmin = (size(block1,2)/2 + size(block3,2)/2)*dt + block2_minDuration;

% Max. area of GY_PE
GY_ramp         = 0;
GY_max          = 0.95*G_max;
ky(1,1)         = ((acquiredkspace(1,2)-1)/2)*Dky;
GY_area(1,1) 	= ky(1,1)/gamma; 

% Max. area of GZ_PE
GZ_PE_ramp          = 0;
GZ_PE_max           = 0.95*G_max;
kz(1,1)             = ((acquiredkspace(1,3)-1)/2)*Dkz;
GZ_PE_area(1,1) 	= kz(1,1)/gamma; 

if minimizeTE == 1
    GZ_PE_duration = block2_minDuration;
else
    GZ_PE_duration = block2_minDuration-size(GZ_ref,2)*dt;
end

while (GY_area(1,1)/block2_minDuration > GY_max) || ...
        (GZ_PE_area(1,1)/GZ_PE_duration > GZ_PE_max)

    block2_minDuration  = block2_minDuration + dt;
    TEmin               = TEmin + dt;
    flag                = 1;
    
    if minimizeTE == 1
        GZ_PE_duration = block2_minDuration;
    else
        GZ_PE_duration = block2_minDuration-size(GZ_ref,2)*dt;
    end

end

if TE < TEmin
    
    if flag == 1
        
        msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or the maximum gradient strength (this option',...
             ' may not be available).'];        
        
    else
    
        msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or decreasing the duration of the RF pulse or the duration of',...
            ' the readout (this option may not be available).'];
        
    end
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end

%% 2nd Block - PE

block2_Duration = TE - (size(block1,2)/2 + size(block3,2)/2)*dt;

% STEP 6.1 - PE
if GY_area(1,1)/block2_Duration <= GY_max
    
    if GY_area(1,1)/GY_max < dt
        
        msg = ['The selected dt is larger than the needed dt (',...
            num2str(GY_area(1,1)/GY_max),' sec) for running the max. PE ',...
            ' gradient. Consider decreasing the dt.'];
        
        error(msg)
        
    end

    GY_magn = GY_area(1,1)/block2_Duration;
    
    [GY,GY_plateau,GY_area_1] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,block2_Duration,0,'GY1',dt);
    

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

% STEP 6.2 - Slice encoding
if GZ_PE_area(1,1)/block2_Duration <= GZ_PE_max
    
    if minimizeTE == 1
        GZ_PE_duration = block2_Duration;
    else
        GZ_PE_duration = block2_Duration-size(GZ_ref,2)*dt;
    end
    
%     if GZ_PE_area(1,1)/GZ_PE_max < dt
%         
%         msg = ['The selected dt is larger than the needed dt (',...
%             num2str(GZ_PE_area(1,1)/GZ_PE_max),' sec) for running the max. Slice Encoding ',...
%             ' gradient. Consider decreasing the dt.'];
%         
%         error(msg)
%         
%     end

    GZ_PE_magn = GZ_PE_area(1,1)/GZ_PE_duration;
    
    [GZ_PE,GZ_PE_plateau,GZ_PE_area_1] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GZ_PE_magn,GZ_PE_duration,0,'GZ_PE1-1',dt);
    

    %------CONTROL-------%
    % Check if the area of the Gy_max is finally equal to the wanted area
    % of Gy_area(1,n). I am using the str2num(num2str()) convertion since
    % isequal() doesn't work
    GZ_PE_area_total = sum(GZ_PE_area_1,2);
    if ~isequal(str2num(num2str(GZ_PE_area(1,1))),...
            str2num(num2str(GZ_PE_area_total))) %#ok
        disp('GZ_PE: PROBLEM WITH THE AREA OF maxGZ_PE');
    end
    %------END CONTROL-------%    
    
end

block2 = zeros(8,floor(block2_Duration/dt));
block2(4,end-size(GX_pre,2)+1:end) = GX_pre;
block2(5,:) = GY;
if minimizeTE==1
    block2_GZ = [GZ_ref+GZ_PE(1,1:size(GZ_ref,2)),GZ_PE(size(GZ_ref,2)+1:end)];
else
    block2_GZ = [GZ_ref,GZ_PE];
end
block2(6,1:size(block2_GZ,2)) = block2_GZ;

%% BUILD TR
TRmin = size(block1,2)+size(block2,2)+size(block3,2);

if N_TR<TRmin
    
    msg = ['The selected TR is lower that the minimum available TR (',...
            num2str(TRmin*dt),' sec) for the current configuration. Consider ',...
            'increasing the TR.'];
        
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end
%% Inversion RF pulse

% Since TI counts from the center of the IR pulse till the center of the
% first excitation pulse, we have to update the time between the end of the
% IR pulse till the start of the first excitation pulse.
TIendPoint          = round(rf_dur/2/dt);

structIR.TI         = str2num(struct_pulseq.ir_time);
structIR.type       = struct_pulseq.ir_type;
structIR.IRduration = str2num(struct_pulseq.ir_duration);
structIR.angle      = str2num(struct_pulseq.preparation_rf_angle);
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
    structExper,TIendPoint,'MP-RAGE',compressIRmodule,conn_localdb,...
    experiment_id,pulseq_id);

%% Synthesize the MPRAGE-TR

pulse_sequence_TR = zeros(8,N_TR);
pulse_sequence_TR(:,1:size(block1,2)) = block1;
pulse_sequence_TR(:,size(block1,2)+1:...
    size(block1,2)+size(block2,2)) = block2;
pulse_sequence_TR(:,size(block1,2)+...
    size(block2,2)+1:size(block1,2)+size(block2,2)+size(block3,2)) = block3;

soft_crushers_TR = zeros(1,N_TR);
soft_crushers_TR(1,N_TR) = 1;

pulse_sequence_MPRAGE_TR    = repmat(pulse_sequence_TR,1,acquiredkspace(1,3));
soft_crushers_MPRAGE_TR     = repmat(soft_crushers_TR,1,acquiredkspace(1,3));

pulse_sequence_MPRAGE_TR(7,:) = 1:size(pulse_sequence_MPRAGE_TR,2);
pulse_sequence_MPRAGE_TR(8,:) = 1:size(pulse_sequence_MPRAGE_TR,2);

% Add the blockIR in front of the pulse_sequence_MPRAGE_TR
pulse_sequence_MPRAGE_TR(7,:) = blockIR_actual_timesteps + (1:size(pulse_sequence_MPRAGE_TR,2));
pulse_sequence_MPRAGE_TR(8,:) = blockIR_actual_timesteps + (1:size(pulse_sequence_MPRAGE_TR,2));
pulse_sequence_MPRAGE_TR      = [blockIR,pulse_sequence_MPRAGE_TR];
% soft_crushers_MPRAGE_TR       = [zeros(1,size(blockIR,2)),soft_crushers_MPRAGE_TR];

% Add one more timepoint at the end of pulse_sequence_MPRAGE_TR that holds
% the total time effect of the TD part of the pulse sequence
pulse_sequence_MPRAGE_TR = [pulse_sequence_MPRAGE_TR,zeros(8,1)];
pulse_sequence_MPRAGE_TR(7,end) = pulse_sequence_MPRAGE_TR(8,end-1)+1;
pulse_sequence_MPRAGE_TR(8,end) = N_MPRAGETR;

% add one extra zero at the end due to the extra point we added on the 
% pulse_sequence_MPRAGE_TR
soft_crushers_MPRAGE_TR = [zeros(1,size(blockIR,2)),soft_crushers_MPRAGE_TR,0];

% Add SwC at the end of blockIR
soft_crushers_MPRAGE_TR(1,size(blockIR,2)) = 1;

pulse_sequence  = repmat(pulse_sequence_MPRAGE_TR,1,acquiredkspace(1,2));
soft_crushers   = repmat(soft_crushers_MPRAGE_TR,1,acquiredkspace(1,2));

isInKspace = zeros(1,size(pulse_sequence,2));

kspace_times    = zeros(acquiredkspace(1,2)*acquiredkspace(1,3),2);
TR_times        = zeros(acquiredkspace(1,2)*acquiredkspace(1,3),2);
MPRAGE_TR_times = zeros(acquiredkspace(1,2),2);

temp1 = 0;
index_kspace = 0;
for i=1:acquiredkspace(1,2)
    for j=1:acquiredkspace(1,3)
        index_kspace = index_kspace + 1;
        kspace_times(index_kspace,1) = ((i-1)*size(pulse_sequence_MPRAGE_TR,2)) + ...
            size(blockIR,2) + (j-1)*N_TR + size(block1,2) + size(block2,2) + 1;
        kspace_times(index_kspace,2) = kspace_times(index_kspace,1) + (GX_plateau(1,3)-GX_plateau(1,2));

        TR_times(index_kspace,1) = ((i-1)*size(pulse_sequence_MPRAGE_TR,2)) + ...
            size(blockIR,2) + (j-1)*N_TR + 1;
        TR_times(index_kspace,2) = ((i-1)*size(pulse_sequence_MPRAGE_TR,2)) + ...
            size(blockIR,2) + j*N_TR;

        isInKspace(1,kspace_times(index_kspace,1):kspace_times(index_kspace,2)) = temp1 + [1:GX_plateau(1,3)];
        temp1 = index_kspace*GX_plateau(1,3);
    end
    
    MPRAGE_TR_times(i,1) = ((i-1)*size(pulse_sequence_MPRAGE_TR,2)) + 1;
    MPRAGE_TR_times(i,2) = i*size(pulse_sequence_MPRAGE_TR,2);
    
    if i==1
        continue;
    else
        pulse_sequence(7,MPRAGE_TR_times(i,1):MPRAGE_TR_times(i,2)) = ...
            pulse_sequence(8,MPRAGE_TR_times(i-1,1)) + ...
            pulse_sequence_MPRAGE_TR(7,:);
        pulse_sequence(8,MPRAGE_TR_times(i,1):MPRAGE_TR_times(i,2)) = ...
            pulse_sequence(8,MPRAGE_TR_times(i-1,1)) + ...
            pulse_sequence_MPRAGE_TR(8,:);
    end
end

kspaceLine = 0;
for n=1:acquiredkspace(1,2)
    for m=1:acquiredkspace(1,3)
        kspaceLine = kspaceLine + 1;       
        
        if ((acquiredkspace(1,2)-1)/2-(n-1))==0   % CONTROL: In case we have an odd 
            % number of kspace lines, the Gy for the line that passes through 
            % zero must have the following characteristics
            ky(1,n) = 0;
            GY_area(1,n) = 0; 
            GY_magn(1,kspaceLine) = 0;
        else
            ky(1,n)=((acquiredkspace(1,2)-1)/2-(n-1))*Dky;
            GY_area(1,n) = ky(1,n)/gamma;  
            GY_magn(1,kspaceLine) = GY_area(1,n)/block2_Duration;
        end
        
        if ((acquiredkspace(1,3)-1)/2-(m-1))==0   % CONTROL: In case we have an odd 
            % number of kspace lines, the Gy for the line that passes through 
            % zero must have the following characteristics
            kz(1,m) = 0;
            GZ_PE_area(1,m) = 0; 
            GZ_PE_magn(1,kspaceLine) = 0;
        else
            kz(1,m)=((acquiredkspace(1,3)-1)/2-(m-1))*Dkz;
            GZ_PE_area(1,m) = kz(1,m)/gamma;  
            GZ_PE_magn(1,kspaceLine) = GZ_PE_area(1,m)/GZ_PE_duration;
        end

        [GY,GY_plateau,GY_area_2] = ...
            pulseSequenceGenerator.addGRtrapezoid(0,...
            GY_magn(1,kspaceLine),block2_Duration,0,['GY',num2str(n),'-',num2str(m)],dt);
        
        [GZ_PE,GZ_PE_plateau,GZ_PE_area_2] = ...
            pulseSequenceGenerator.addGRtrapezoid(0,...
            GZ_PE_magn(1,kspaceLine),GZ_PE_duration,0,['GZ_PE',num2str(n),'-',num2str(m)],dt);

        %------CONTROL-------%
        % Check if the area of the GY_max and GZ_PE_max are  equal to the 
        % wanted area of GY_area(1,n) and GZ_PE_area(1,m). I am using the 
        % str2num(num2str()) convertion since isequal() doesn't work
        GY_area_total = sum(GY_area_2,2);
        if ~isequal(str2num(num2str(GY_area(1,n))),...
                str2num(num2str(GY_area_total))) %#ok
            disp(['Gy: PROBLEM WITH THE AREA OF GY',num2str(n)]);
        end
        GZ_PE_area_total = sum(GZ_PE_area_2,2);
        if ~isequal(str2num(num2str(GZ_PE_area(1,m))),...
                str2num(num2str(GZ_PE_area_total))) %#ok
            disp(['Gz: PROBLEM WITH THE AREA OF GZ_PE',num2str(m)]);
        end
        %------END CONTROL-------%

        pulse_sequence(5,TR_times(kspaceLine,1)+size(block1,2)+1:...
            TR_times(kspaceLine,1)+size(block1,2)+size(GY,2)) = GY;
        
        if minimizeTE==1
            block2_GZ = [GZ_ref+GZ_PE(1,1:size(GZ_ref,2)),GZ_PE(size(GZ_ref,2)+1:end)];
        else
            block2_GZ = [GZ_ref,GZ_PE];
        end
        pulse_sequence(6,TR_times(kspaceLine,1)+size(block1,2)+1:...
            TR_times(kspaceLine,1)+size(block1,2)+size(block2_GZ,2)) = block2_GZ;        
    end
end

%% Calculate extra parameters
encodedkSpace.lines         = acquiredLines;
encodedkSpace.columns       = round(GX_plateau_duration/dt);
encodedkSpace.contrasts     = 1;  % since we acquire only one echoe per TR

N_pulse = size(pulse_sequence,2);

PulseqDuration = MPRAGETR*acquiredkspace(1,2);
disp(['The total duration of the pulse sequence is ',num2str(PulseqDuration),'sec'])

%% FINAL CHECK
if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
        max(pulse_sequence(6,:))>=G_max
    exception = MException('GrMagn:GreaterThanMax',...
        'The gradients magnitude exceeds the maximum gradient strength');
    throw(exception)
end