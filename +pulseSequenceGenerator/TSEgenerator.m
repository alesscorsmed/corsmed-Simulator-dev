function [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,...
    N_TR,kspace_times,BW1,TR_times,TRmin_tmstps,orderOfKspacelines,...
    TEeffective,PulseqDuration,encodedkSpace] = TSEgenerator(RFmatrix1,...
    RFmatrix2,acquiredkspace,FOV,receiver_BW,G_max,TR,TE,dt,gamma,ETL,...
    conn_localdb,experiment_id,pulseq_id,struct_pulseq,pulseSeqFamilyName,...
    partialFourierStruct)
% This function generates Turbo Spin Echo (TSE) pulse sequences. It does 
% NOT take into account GR ramps. ETL stands for Echo Train Length.

N_TR                = floor(TR/dt);   % Number of time-steps during the TR

rf_dur1             = RFmatrix1(1,1);
cycles1             = RFmatrix1(2,1);
angle1              = RFmatrix1(3,1);
slice_thickness1    = RFmatrix1(4,1);

Dkx                 = 1/FOV(1,1);   %k-space Dkx   ALLAXE TO
Dky                 = 1/FOV(1,2);   %k-space Dky   ALLAXE TO

% Acquire less k-space lines if halfFourierFactor < 1
partialFourierFactor    = partialFourierStruct.factor;
acquiredLines           = round(partialFourierFactor*acquiredkspace(1,2));

% The phaseConjugate should only apply for SS-FSE
if strcmp(partialFourierStruct.type,'phaseConjugate')
    if strcmp(pulseSeqFamilyName,'SS-FSE')
        ETL_original    = ETL;
        ETL             = acquiredLines;
        disp(['ETL changed from ',num2str(ETL_original),' to ',num2str(ETL)])
    end
    % The acquired lines should be an integer multiple of ETL
    multfactor      = ceil(acquiredLines/ETL);
    acquiredLines   = multfactor*ETL;
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
block1(2,GZ_plateau1(1,2):GZ_plateau1(1,3)) = 90*pi/180; % phase 90 (in rads)

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

disp(['Gss 90RF (area)      : ',num2str(GZ_area1)])
disp(['Gss 180RF (area)     : ',num2str(GZ_area2)])

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

%% 2nd Block - GR_Z (90) ref. + GR_Z (180) ref. + GR_X pre.

% STEP 3 - GZ refocussing gradient (GZ_ref) and crushers on either sides of
% the Gz for the 180 RF pulse
% We assume that the refocusing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "slice selection" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GZ_area(1,2)/2 + GZ_area(1,3))
GZ2_ref_magn             = 0.95*G_max;                                                                                    
Gz2_ref_plateau_tmstps   = floor((GZ_area2(1,2)/2 + GZ_area2(1,3))/...
    (GZ2_ref_magn*dt));
G2Z_ref_magn             = (GZ_area2(1,2)/2 + GZ_area2(1,3))/...
    (Gz2_ref_plateau_tmstps*dt);  % recalculate GZ_ref_magn so as to yield 
        % equal gradient areas for the time Gz_ref_plateau_tmstps

[GZ2_ref,GZ2_ref_plateau,GZ2_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    GZ2_ref_magn,Gz2_ref_plateau_tmstps*dt,0,'GZ2_ref',dt);

disp(['GssRef 180RF (area)  : ',num2str(GZ2_ref_area)])

[GZ_ref,GZ_ref_plateau,GZ_ref_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,...
    GZ2_ref_magn,Gz2_ref_plateau_tmstps*dt,0,'GZ_ref',dt);
GZ_ref = -GZ_ref;
GZ_ref_area = -GZ_ref_area;

disp(['GssRef 90RF (area)   : ',num2str(GZ_ref_area)])

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
GX_pre_plateau_tmstps = floor((GxMainPart + GX_area(1,1))/(GX_pre_magn*dt));
GX_pre_magn = (GxMainPart + GX_area(1,1))/(GX_pre_plateau_tmstps*dt);

[GX_pre,GX_pre_plateau,GX_pre_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_pre_magn,...
    GX_pre_plateau_tmstps*dt,0,'GX_pre',dt);

%% Calculate the min TE
block2_minDuration = (max(size(GZ_ref,2)+size(GZ2_ref,2),size(GX_pre,2)))*dt;
TEmin_tmstps = 2*((size(block1,2)/2 + size(block3,2)/2)) + ...
    max(size(GZ_ref,2)+size(GZ2_ref,2),size(GX_pre,2));
TEmin = 2*((size(block1,2)/2 + size(block3,2)/2)*dt + block2_minDuration);

% If there is an overlap between the second half of the second RF and the
% first (till TE) part of readout
if TEmin/2 < (size(block3,2)/2 + size(GZ2_ref,2) + (1 - 0.5/readConjugateFactor)*size(GX,2))*dt
    TEminProposed = 2*(size(block3,2)/2 + size(GZ2_ref,2) + (1 - 0.5/readConjugateFactor)*size(GX,2))*dt;
    msg = ['The selected TE is smaller than the minimum available TE (',...
        num2str(TEminProposed),'sec) for the current configuration. Consider',...
        ' increasing the TE or decreasing the duration of the readout.'];
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end

% Available space for GY_PE
availGYspace_tmstps = TEmin_tmstps - ...
    (size(block1,2)/2 + size(block3,2) + max(size(GZ_ref,2)+size(GZ2_ref,2),size(GX_pre,2))) - ...
    size(GZ2_ref,2) - round((1 - 0.5/readConjugateFactor)*size(GX,2));  % @@@ Make a change here

while availGYspace_tmstps <= 0
    
    TEmin_tmstps        = TEmin_tmstps + 1;
    TEmin               = TEmin + dt;
    
    availGYspace_tmstps = TEmin_tmstps - ...
        (size(block1,2)/2 + size(block3,2) + max(size(GZ_ref,2)+size(GZ2_ref,2),size(GX_pre,2))) - ...
        size(GZ2_ref,2) - round((readConjugateFactor - 0.5)*size(GX,2)); % @@@ Make a change here
    
end

availGYspace = availGYspace_tmstps*dt;

% Max. area of GY_PE
GY_ramp         = 0;
GY_max          = 0.95*G_max;
ky(1,1)         = ((acquiredkspace(1,2)-1)/2)*Dky;
GY_area(1,1) 	= ky(1,1)/gamma; 

while GY_area(1,1)/availGYspace > GY_max
    
    TEmin_tmstps        = TEmin_tmstps + 1;
    TEmin               = TEmin + dt;
    availGYspace_tmstps = availGYspace_tmstps + 1;
    availGYspace        = availGYspace_tmstps*dt;
    flag                = 1;

end

if TE < TEmin
    
    if flag == 1
        
        msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or the maximum gradient strength (if applicable).'];        
        
    else
    
        msg = ['The selected TE is smaller than the minimum available TE (',...
            num2str(TEmin),'sec) for the current configuration. Consider',...
            ' increasing the TE or decreasing the duration of the RF pulse or the duration of',...
            ' the readout (this option may not be available).'];
        
    end
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end

%% 2nd Block - Build this block
% TE is equal to 2*(block1/2+block2+block3/2) => 
% block2 = (TE - block1 - block3)/2
block2_Duration = (TE - (size(block1,2) + size(block3,2))*dt)/2;
block2          = zeros(8,floor(block2_Duration/dt));

block2(6,1:size(GZ_ref,2))          = GZ_ref;
block2(4,end-size(GX_pre,2)+1:end)  = GX_pre;
block2(6,end-size(GZ2_ref,2)+1:end)  = GZ2_ref;

%% 4th Block - Add PE and build the block

% STEP 6 - PE
if GY_area(1,1)/availGYspace <= GY_max
    
    if GY_area(1,1)/GY_max < dt
        
        msg = ['The selected dt is larger than the needed dt (',...
            num2str(GY_area(1,1)/GY_max),' sec) for running the max. PE ',...
            ' gradient. Consider decreasing the dt.'];
        
        error(msg)
        
    end

    GY_magn = GY_area(1,1)/availGYspace;
    
    [GY,GY_plateau,GY_area_1] = ...
        pulseSequenceGenerator.addGRtrapezoid(GY_ramp,...
        GY_magn,availGYspace,GY_ramp,'GY1',dt);
    

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

% According to figure 16.40, pg.777 of the Handbook of pulse sequences:
% block1 + block2 + block3 + block4 = 3*tau = 3*(TE/2) => 
% block4 = 1.5*TE - block1 - block2 - block3
block4_duration_tmstps  = 1.5*round(TE/dt)-size(block1,2)-...
    size(block2,2)-size(block3,2);
block4                  = zeros(8,block4_duration_tmstps);

current_tmstp_start = 1;
current_tmstp_end   = size(GZ2_ref,2);
block4(6,current_tmstp_start:current_tmstp_end) = GZ2_ref;

current_tmstp_start = size(GZ2_ref,2)+1;
current_tmstp_end   = size(GZ2_ref,2)+size(GY,2);
block4(5,current_tmstp_start:current_tmstp_end) = GY;

current_tmstp_start = round(block4_duration_tmstps/2)-round(size(GX,2)/2)+1;
current_tmstp_end   = round(block4_duration_tmstps/2)-round(size(GX,2)/2)+size(GX,2);
block4(4,current_tmstp_start:current_tmstp_end) = GX;

current_tmstp_end   = block4_duration_tmstps;
current_tmstp_start = block4_duration_tmstps - size(GZ2_ref,2) + 1;
block4(6,current_tmstp_start:current_tmstp_end) = GZ2_ref;

current_tmstp_end   = current_tmstp_start - 1;
current_tmstp_start = current_tmstp_end - size(GY,2) + 1;
block4(5,current_tmstp_start:current_tmstp_end) = -GY;

%% BUILD TR
TRmin_tmstps = size(block1,2) + size(block2,2) + ...
    ETL*(size(block3,2) + size(block4,2));

if N_TR<TRmin_tmstps
    
    msg = ['The selected TR is lower than the minimum available TR (',...
            num2str(TRmin_tmstps*dt),' sec) for the current configuration. Consider ',...
            'increasing the TR or decreasing the echo train length (ETL).'];
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end

noOfTRs = ceil(acquiredLines/ETL);

pulse_sequence              = [];
soft_crushers               = [];

TRendTmstp                  = 0;

kspace_times                = zeros(acquiredLines,2);
TR_times                    = zeros(acquiredLines,2);
GY_times                    = zeros(acquiredLines,2);
GY_times_refoc              = zeros(acquiredLines,2);

for iTR = 1:noOfTRs
    
    % If fatsat in on
    if strcmp(pulseSeqFamilyName,'IR-TSE') || ...
            ((strcmp(pulseSeqFamilyName,'TSE') || strcmp(pulseSeqFamilyName,'SS-FSE')) && ...
            isfield(struct_pulseq,'fatsat') && ...
            strcmp(struct_pulseq.fatsat,'lipidir'))

        structIR.TI         = str2num(struct_pulseq.ir_time);

        structIR.angle      = 180;
        structIR.phase      = 0; % in degrees
        structIR.IRcycles   = 2; % only for sinc RF pulse

        structExper.dt      = dt;
        structExper.gamma   = gamma;

        compressIRmodule    = 1;

        if ((strcmp(pulseSeqFamilyName,'TSE') || strcmp(pulseSeqFamilyName,'SS-FSE')) && ...
                isfield(struct_pulseq,'fatsat') && ...
                strcmp(struct_pulseq.fatsat,'lipidir'))
            structIR.type       = 'sinc';
            structIR.IRduration = 0.003;
        elseif strcmp(pulseSeqFamilyName,'IR-TSE')
            structIR.type       = struct_pulseq.ir_type;
            structIR.IRduration = str2num(struct_pulseq.ir_duration);
        end

        TIendPoint = round(size(RF1,2)/2);

        if iTR == 1 % run this the very first time

            % blockIR_actual_timesteps holds the "true" number of timesteps of the
            % blockIR. blockIR_timesteps holds the "effective" number of timesteps due
            % to the application of the fast algorithm

            [blockIR,blockIR_actual_timesteps,blockIR_timesteps] = ...
                pulseSequenceGenerator.generateIRblock(structIR,...
                structExper,TIendPoint,'TSE',compressIRmodule,conn_localdb,...
                experiment_id,pulseq_id);
        end

        blockIR(7,:)  = TRendTmstp+blockIR(7,:);
        blockIR(8,:)  = TRendTmstp+blockIR(8,:);

        pulse_sequence          = [pulse_sequence,blockIR];
        soft_crushers           = [soft_crushers,zeros(1,size(blockIR,2))];

        TRendTmstp              = TRendTmstp + blockIR_timesteps;
    end
    
    pulse_sequence_TR = zeros(8,TRmin_tmstps);
    pulse_sequence_TR(:,1:size(block1,2)) = block1;
    pulse_sequence_TR(:,size(block1,2)+1:...
        size(block1,2)+size(block2,2)) = block2;
    current_tmstp = size(block1,2) + size(block2,2);

    for iEcho = 1:ETL
        pulse_sequence_TR(:,current_tmstp+1:...
            current_tmstp+size(block3,2)+size(block4,2)) = [block3,block4];        
        
        TR_times((iTR-1)*ETL+iEcho,1) = TRendTmstp + current_tmstp + 1;
        TR_times((iTR-1)*ETL+iEcho,2) = TRendTmstp + current_tmstp + size(block3,2) + size(block4,2);
        GY_times((iTR-1)*ETL+iEcho,1) = TRendTmstp + current_tmstp + size(block3,2) + size(GZ2_ref,2) + 1;
        GY_times((iTR-1)*ETL+iEcho,2) = TRendTmstp + current_tmstp + size(block3,2) + size(GZ2_ref,2) + size(GY,2);
        kspace_times((iTR-1)*ETL+iEcho,1) = TRendTmstp + current_tmstp + size(block3,2) + round(block4_duration_tmstps/2) - round(size(GX,2)/2) + 1;
        kspace_times((iTR-1)*ETL+iEcho,2) = TRendTmstp + current_tmstp + size(block3,2) + round(block4_duration_tmstps/2) - round(size(GX,2)/2) + size(GX,2);
        GY_times_refoc((iTR-1)*ETL+iEcho,1) = TRendTmstp + current_tmstp + size(block3,2) + size(block4,2) - size(GZ2_ref,2) - size(GY,2) + 1;
        GY_times_refoc((iTR-1)*ETL+iEcho,2) = TRendTmstp + current_tmstp + size(block3,2) + size(block4,2) - size(GZ2_ref,2);
        
        current_tmstp = current_tmstp+size(block3,2)+size(block4,2);
    end

    pulse_sequence_TR(7,:) = TRendTmstp+[1:size(pulse_sequence_TR,2)];
    pulse_sequence_TR(8,:) = TRendTmstp+[1:size(pulse_sequence_TR,2)];
    
    % Add two more extra points. The first holds the accumulated effect of the
    % pulse sequence till the end of the TR (minus one dt) and the second holds
    % the temporal effect of the soft crushers at the end of the TR
    pulse_sequence_TR           = [pulse_sequence_TR,[zeros(6,1);...
        TRendTmstp+TRmin_tmstps+1;TRendTmstp+N_TR-1]];
    pulse_sequence_TR           = [pulse_sequence_TR,[zeros(6,1);...
        TRendTmstp+N_TR;TRendTmstp+N_TR]];

    soft_crushers_TR            = zeros(1,size(pulse_sequence_TR,2));
    soft_crushers_TR(1,size(pulse_sequence_TR,2)) = 1;
    
    pulse_sequence              = [pulse_sequence,pulse_sequence_TR];
    soft_crushers               = [soft_crushers,soft_crushers_TR];
    
    TRendTmstp = TRendTmstp + size(pulse_sequence_TR,2);
    
end

isInKspace = zeros(1,size(pulse_sequence,2));

temp1 = 0;
for i=1:acquiredLines   
    isInKspace(1,kspace_times(i,1):kspace_times(i,2)) = temp1 + [1:GX_plateau(1,3)];
    temp1 = i*GX_plateau(1,3);
end

% Define the order of acquisition of the kspace lines
orderOfKspace = 'segmental';
if strcmp(orderOfKspace,'sequential')
    orderOfKspacelines = 1:acquiredLines;
elseif strcmp(orderOfKspace,'segmental')
    orderOfKspacelines = zeros(1,acquiredLines);
    linesPerSegment = acquiredLines/ETL;
    currentSegment = 1;
    currentLineWithinSegment = 1;
    for i=1:acquiredLines
        orderOfKspacelines(1,i) = (currentSegment-1)*linesPerSegment+currentLineWithinSegment;
        if currentSegment<ETL
            currentSegment = currentSegment + 1;
        else
            currentSegment = 1;
            currentLineWithinSegment = currentLineWithinSegment + 1;
        end
    end
end

iechoe = 1;
for n = 2:size(orderOfKspacelines,2) % kspace(1,2)
    
    iechoe = iechoe + 1;
    kspace_line = orderOfKspacelines(1,n);
    
    if ((acquiredkspace(1,2)-1)/2-(kspace_line-1))==0   % CONTROL: In case we have an odd 
        % number of kspace lines, the Gy for the line that passes through 
        % zero must have the following characteristics
        ky(1,kspace_line) = 0;
        GY_area(1,kspace_line) = 0; 
        GY_magn = 0;
    else
        ky(1,kspace_line)=((acquiredkspace(1,2)-1)/2-(kspace_line-1))*Dky;
        GY_area(1,kspace_line) = ky(1,kspace_line)/gamma;  
        GY_magn = GY_area(1,kspace_line)/availGYspace;
    end
    
    [GY,GY_plateau,GY_area_2] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,availGYspace,0,['GY',num2str(kspace_line)],dt);
                                           
    %------CONTROL-------%
    % Check if the area of the Gy_max is finally equal to the wanted area
    % of Gy_area(1,n). I am using the str2num(num2str()) convertion since
    % isequal() doesn't work
    GY_area_total = sum(GY_area_2,2);
    if ~isequal(str2num(num2str(GY_area(1,kspace_line))),...
            str2num(num2str(GY_area_total))) %#ok
        disp(['Gy: PROBLEM WITH THE AREA OF GY',num2str(kspace_line)]);
    end
    %------END CONTROL-------%
    
    pulse_sequence(5,GY_times(iechoe,1):GY_times(iechoe,2)) = GY;
    pulse_sequence(5,GY_times_refoc(iechoe,1):GY_times_refoc(iechoe,2)) = -GY;
    
end

N_pulse = size(pulse_sequence,2);

%% Calculate extra parameters
encodedkSpace.lines         = acquiredLines;
encodedkSpace.columns       = round(GX_plateau_duration/dt);

if strcmp(pulseSeqFamilyName,'SS-FSE')
    TEeffective = (acquiredkspace(1,2)/2)*TE+TE/2;
else
    TEeffective = (ETL/2)*TE+TE/2;
end
disp(['The effective TE is ',num2str(TEeffective),'sec'])

if exist('N_IRmodule','var')
    PulseqDuration = noOfTRs*(TR+blockIR_actual_timesteps*dt);
else
    if isfield(struct_pulseq,'fatsat') && strcmp(struct_pulseq.fatsat,'optimal')
        if isfield(struct_pulseq,'ir_duration')
            IRduration = str2num(struct_pulseq.ir_duration);
        else
            IRduration = 0.003;
        end
        PulseqDuration = noOfTRs*(TR+str2num(struct_pulseq.ir_time)+...
            IRduration/2-rf_dur1/2);
    else
        PulseqDuration = noOfTRs*TR;
    end
end
disp(['The total duration of the pulse sequence is ',num2str(PulseqDuration),'sec'])

%% FINAL CHECK
if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
        max(pulse_sequence(6,:))>=G_max
    exception = MException('GrMagn:GreaterThanMax',...
        ['The gradients magnitude (X:',num2str(max(pulse_sequence(4,:))),...
        ', Y:',num2str(max(pulse_sequence(5,:))),') exceeds the maximum',...
        ' gradient strength (',num2str(G_max),')']);
    throw(exception)
end

%% PLOT
% timeaxis = [1:N_pulse]*dt;
% 
% a1 = subplot(7,1,1);
% plot(timeaxis,pulse_sequence(1,:),'r')
% ylabel('RF magnitude') % in Tesla
% 
% a2 = subplot(7,1,2);
% plot(timeaxis,pulse_sequence(2,:),'r')
% ylabel('RF phase') % in radians
% 
% a3 = subplot(7,1,3);
% plot(timeaxis,pulse_sequence(3,:),'r')
% ylabel('RF Freq.') % in Hz
% 
% a4 = subplot(7,1,4);
% plot(timeaxis,pulse_sequence(4,:))
% ylabel('Gx') % in T/m
% 
% a5 = subplot(7,1,5);
% plot(timeaxis,pulse_sequence(5,:))
% ylabel('Gy') % in T/m
% hold on
% plot(GY_times(:,1)*dt,pulse_sequence(5,GY_times(:,1)'),'rs')
% 
% a6 = subplot(7,1,6);
% plot(timeaxis,pulse_sequence(6,:))
% ylabel('Gz') % in T/m
% 
% a7 = subplot(7,1,7);
% plot(timeaxis,isInKspace,'k')
% ylabel('isInKspace')
% set(gcf,'Name','Pulse Sequence')
% 
% linkaxes([a1 a2 a3 a4 a5 a6 a7],'x');
% set(gca,'XTick',[])