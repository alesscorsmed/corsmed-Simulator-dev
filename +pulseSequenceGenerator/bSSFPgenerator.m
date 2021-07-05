function [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,...
    N_TR,kspace_times,BW,TR_times,receiver_phase,TRmin,PulseqDuration,...
    bSSFPCenterOfAcquisition,encodedkSpace] = ...
    bSSFPgenerator(RFmatrix,acquiredkspace,FOV,receiver_BW,G_max,TR,TE,dt,gamma,...
    bSSFPpreparation,experiment_id,pulseq_id,conn_localdb,partialFourierStruct,...
    lengthRamp,addTailAby2)
% This function generates bSSFP pulse sequences. It does NOT take into
% account ramps.

if TE~=TR/2
    msg = ['The echo time (TE) should be half the repetition time (TR) in ',...
        'the bSSFP pulse sequence design.'];

    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end

N_TR            = floor(TR/dt);   % Number of time-steps during the TR

rf_dur          = RFmatrix(1,1);
cycles          = RFmatrix(2,1);
angle           = RFmatrix(3,1);
slice_thickness = RFmatrix(4,1);

Dkx             = 1/FOV(1,1);   %k-space Dkx   ALLAXE TO
Dky             = 1/FOV(1,2);   %k-space Dky   ALLAXE TO

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

N_pulse         = acquiredLines*N_TR;
G_pulse         = zeros(3,N_pulse);
B1_pulse        = zeros(2,N_pulse);

flag            = 0;

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

% The phase of the first RF pulse is -90. The same goes for the receiver.
block1(2,GZ_plateau(1,2):GZ_plateau(1,3)) = -90*pi/180; % phase -90 (in rads)

receiver_phase(1,1) = -90*pi/180; % phase -90 (in rads)


%% 3rd Block - GR_X readout

% STEP 5 - GR_X
GX_magn                 = receiver_BW/(gamma*FOV(1,1));
GX_plateau_duration     = acquiredRO/receiver_BW;  
[GX,GX_plateau,GX_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_magn,...
    GX_plateau_duration,0,'GX',dt);

block3 = zeros(8,size(GX,2));
block3(4,1:size(GX,2)) = GX;

%% 2nd Block - GR_Z ref. + GR_X pref.

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

% STEP 4.1 - GR_X pre.
% We assume that the pre-phasing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,1)) if no partial echo is selected.
% Otherwise is calculated taking into consideration the area before the TE.

% GX_pre_magn = 0.95*G_max;
GX_pre_ref_tmstps = ceil(sum(GX)/(0.90*G_max));
GX_pre_magn = sum(GX)/(GX_pre_ref_tmstps);
if strcmp(partialFourierStruct.type,'readConjugate')
    GxMainPart = (readConjugateFactor - 0.5) * (GX_area(1,2)/readConjugateFactor);
    
    
    GX_pre_plateau_tmstps = ceil((readConjugateFactor - 0.5) * GX_pre_ref_tmstps/readConjugateFactor); % SBT210120 fix to ensure balanced gradients

else
    GxMainPart = GX_area(1,2)/2;
    
    
    GX_pre_plateau_tmstps = ceil(GX_pre_ref_tmstps/2); % SBT210120 fix to ensure balanced gradients

end
% GX_pre_plateau_tmstps = ceil((GxMainPart + GX_area(1,1))/(GX_pre_magn*dt)); % SBT210120 changed floor to ceil
% GX_pre_magn = (GxMainPart + GX_area(1,1))/(GX_pre_plateau_tmstps*dt);
GX_pre_magn = GX_pre_magn*(-1);

GX_pre = nan(1, GX_pre_plateau_tmstps);
GX_pre(1,:) = GX_pre_magn;

% [GX_pre,GX_pre_plateau,GX_pre_area] = ...
%     pulseSequenceGenerator.addGRtrapezoid(0,GX_pre_magn,...
%     GX_pre_plateau_tmstps*dt,0,'GX_pre',dt);

% STEP 4.2 - GR_X ref.
% We assume that the refocussing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,3)) if no partial echo is selected.
% Otherwise is calculated taking into consideration the area after the TE.

% GX_ref_magn = 0.95*G_max;
if strcmp(partialFourierStruct.type,'readConjugate')
    GxMainPart2 = 0.5 * (GX_area(1,2)/readConjugateFactor);
else
    GxMainPart2 = GX_area(1,2)/2;
end

GX_ref_plateau_tmstps = GX_pre_ref_tmstps - GX_pre_plateau_tmstps; % SBT210120 fix to ensure balanced gradients
% GX_ref_plateau_tmstps = ceil(sum(GX)/(GX_ref_magn)); % SBT210120 fix to ensure balanced gradients
GX_ref_magn = GX_pre_magn;

GX_ref = nan(1, GX_ref_plateau_tmstps);
GX_ref(1,:) = GX_ref_magn;

% [GX_ref,GX_ref_plateau,GX_ref_area] = ...
%     pulseSequenceGenerator.addGRtrapezoid(0,GX_ref_magn,...
%     GX_ref_plateau_tmstps*dt,0,'GX_ref',dt);

%% 200121 Check read gradients and force perfect balancing:
gx_0thmoment = sum([GX_pre, GX, GX_ref]);
% if (gx_0thmoment > 1e-10)
%     % Unbalanced read-gradients found: Force symmetric pre and rephasors
%     gx_pre_ref_timesteps = numel(GX_pre) + numel(GX_ref);
%     gx_pre_ref_amplitude = sum(GX)/gx_pre_ref_timesteps;
%     GX_pre(1:end) = -gx_pre_ref_amplitude;
%     GX_ref(1:end) = -gx_pre_ref_amplitude;
%     % gx_0thmoment = sum([GX_pre, GX, GX_ref]);
% end

%% Calculate the min TE

block2_minDuration = max([size(GZ_ref,2), size(GX_pre,2)])*dt; % 200121SBT Added missing [] within max() call
% TEmin = (size(block1,2)/2 + ceil((1 - 0.5/readConjugateFactor)*size(block3,2)))*dt + block2_minDuration; % 200121SBT replace round with ceil

% Max. area of GY_PE
GY_ramp         = 0;
GY_max          = 0.95*G_max;
ky(1,1)         = ((acquiredkspace(1,2)-1)/2)*Dky;
GY_area(1,1) 	= ky(1,1)/gamma; 

block2_minDuration_temp = dt*ceil(GY_area(1,1)/(dt*GY_max)); % 200121SBT block2minDur calc without having to use while loop
block2_minDuration = max([block2_minDuration, block2_minDuration_temp]); % 200121SBT block2minDur calc without having to use while loop
TEmin = (size(block1,2)/2 + ceil((1 - 0.5/readConjugateFactor)*size(block3,2)))*dt + block2_minDuration; % 200121SBT moved TEmin calc after final block2_minDuration calc

% while GY_area(1,1)/block2_minDuration > GY_max
%     
%     block2_minDuration  = block2_minDuration + dt;
%     TEmin               = TEmin + dt;
%     flag                = 1;
%     
% end

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

%% 2nd Block - PE

block2_Duration = TE - (size(block1,2)/2 + round((1 - 0.5/readConjugateFactor)*size(block3,2)))*dt;
block4_Duration = TE - (size(block1,2)/2 + round((0.5/readConjugateFactor)*size(block3,2)))*dt; 

% STEP 6 - PE
if GY_area(1,1)/block4_Duration <= GY_max
    
    if GY_area(1,1)/GY_max < dt
        
        msg = ['The selected dt is larger than the needed dt (',...
            num2str(GY_area(1,1)/GY_max),' sec) for running the max. PE ',...
            ' gradient. Consider decreasing the dt.'];
        
        error(msg)
        
    end
    
    GY_magn = GY_area(1,1)/block4_Duration;
    
    [GY,GY_plateau,GY_area_1] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,block4_Duration,0,'GY1',dt);
    

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
block2(5,end-size(GY,2)+1:end) = GY;

block4 = zeros(8,floor(block4_Duration/dt));

% Verify that the sum size of the blocks is equal to N_TR
TRmin = size(block1,2) + size(block2,2) + size(block3,2) + size(block4,2);

if N_TR<TRmin
    
    msg = ['The selected TR is lower than the minimum available TR (',...
            num2str(TRmin*dt),' sec) for the current configuration. Consider ',...
            'increasing the TR.'];
        
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
elseif N_TR>TRmin
    
    disp(['WARNING FOR THE ADMINISTRATOR: the total size of the blocks in TR ',...
        'is lower than the selected TR.'])
    
    block4 = [zeros(8,N_TR-TRmin),block4];
    
end

%% 4th block - GRY refoc. + PE refoc. + GRZ refoc.
block4(4,end-size(GX_ref,2)+1:end)  = GX_ref;
block4(5,1:size(GY,2))              = -GY;
block4(6,end-size(GZ_ref,2)+1:end)  = GZ_ref;

%% BUILD TR
pulse_sequence_TR = zeros(8,N_TR);
pulse_sequence_TR(:,1:size(block1,2)) = block1;
pulse_sequence_TR(:,size(block1,2)+1:...
    size(block1,2)+size(block2,2)) = block2;
pulse_sequence_TR(:,size(block1,2)+...
    size(block2,2)+1:size(block1,2)+size(block2,2)+size(block3,2)) = block3;
pulse_sequence_TR(:,size(block1,2)+...
    size(block2,2)+size(block3,2)+1:size(block1,2)+size(block2,2)+...
    size(block3,2)+size(block4,2)) = block4;

soft_crushers_TR = zeros(1,N_TR);
% soft_crushers_TR(1,N_TR) = 1;

pulse_sequence  = repmat(pulse_sequence_TR,1,acquiredLines);
soft_crushers   = repmat(soft_crushers_TR,1,acquiredLines);

isInKspace = zeros(1,size(pulse_sequence,2));

pulse_sequence(5,N_TR+1:end) = ...
    zeros(1,size(N_TR+1:size(pulse_sequence,2),2));

kspace_times    = zeros(acquiredLines,2);
TR_times        = zeros(acquiredLines,2);

temp1 = 0;
for i=1:acquiredLines
    kspace_times(i,1) = ((i-1)*N_TR) + size(block1,2) + size(block2,2) + 1;
    kspace_times(i,2) = kspace_times(i,1) + (GX_plateau(1,3)-GX_plateau(1,2));
    
    TR_times(i,1) = (i-1)*N_TR+1;
    TR_times(i,2) = i*N_TR;
    
    isInKspace(1,kspace_times(i,1):kspace_times(i,2)) = temp1 + [1:GX_plateau(1,3)];
    temp1 = i*GX_plateau(1,3);
end

for n=2:acquiredLines
    if ((acquiredkspace(1,2)-1)/2-(n-1))==0   % CONTROL: In case we have an odd 
        % number of kspace lines, the Gy for the line that passes through 
        % zero must have the following characteristics
        ky(1,n) = 0;
        GY_area(1,n) = 0; 
        GY_magn = 0;
    else
        ky(1,n)=((acquiredkspace(1,2)-1)/2-(n-1))*Dky;
        GY_area(1,n) = ky(1,n)/gamma;  
        GY_magn = GY_area(1,n)/block4_Duration;
    end
    
    [GY,GY_plateau,GY_area_2] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,block4_Duration,0,['GY',num2str(n)],dt);
                                           
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
    
    pulse_sequence(5,((n-1)*N_TR)+size(block1,2)+1:...
        ((n-1)*N_TR)+size(block1,2)+size(GY,2)) = GY;
    
    pulse_sequence(5,((n-1)*N_TR)+size(block1,2)+size(block2,2)+size(block3,2)+1:...
        ((n-1)*N_TR)+size(block1,2)+size(block2,2)+size(block3,2)+size(GY,2)) = -GY;
    
    % Change also the polarity of the RF phase (for every even TR)
    if mod(n,2)==0
        pulse_sequence(2,(n-1)*N_TR+1:n*N_TR) = ...
            -pulse_sequence(2,(n-1)*N_TR+1:n*N_TR);       
    end
    
    sign = cosd(n*180);
    receiver_phase(1,n) = sign*(90*pi/180); % phase (in rads)
    
end


%% ADD PREPULSE AND UPDATE MATRICES
% 5th Block - RF, Gz (slice selection) and Gz refoc.

if ~isempty(bSSFPpreparation)
    if strcmp(bSSFPpreparation,'a/2')

        N_block5 = floor(TR/2/dt);

        % STEP 1 - Design RF
        [RF2,BW2,rf_timesteps]  = pulseSequenceGenerator.addRFsinc(rf_dur,cycles,...
            angle/2,dt,gamma);

        % STEP 2 - Design slice selection gradient (GZ)
        Gz_magn_prepulse        = BW2/(gamma*slice_thickness);
        [GZ_prepulse,GZ_prepulse_plateau,GZ_prepulse_area] = ...
            pulseSequenceGenerator.addGRtrapezoid(0,Gz_magn_prepulse,rf_dur,...
            0,'GZ_prepulse',dt);

        % Check Error: if the length of RF pulse is equal to the length of the plateau
        size_plateau = size(GZ_prepulse_plateau(1,2):GZ_prepulse_plateau(1,3),2);
        if size_plateau~=rf_timesteps
            msg = ['ERROR - DO NOT IGNORE IT - The length of RF pulse is NOT',...
                ' equal to the length of the plateau'];
            disp(msg)
        end

        % STEP 3 - GZ refocussing gradient (GZ_prepulse_ref)
        % We assume that the refocusing pulse is a hard pulse with no ramps. The
        % plateau duration is equal to the half duration of the "slice selection" 
        % pulse and its magnitude is calculated so as the product to be equal to 
        % the (GZ_area(1,2)/2 + GZ_area(1,3))                                                                                     
        Gz_prepulse_ref_plateau_tmstps = N_block5 - size(GZ_prepulse,2);

        % recalculate GZ_prepulse_ref_magn so as to yield equal gradient areas 
        % for the time Gz_prepulse_ref_plateau_tmstps
        GZ_prepulse_ref_magn    = (sum(GZ_prepulse_area,2)/...
            (Gz_prepulse_ref_plateau_tmstps*dt)); 

        if GZ_prepulse_ref_magn>G_max

            if Gz_prepulse_ref_plateau_tmstps <= 0

                msg = ['The duration of the prepulse RF (a/2) is equal or ',...
                    'larger than the TR/2. Consider decreasing the duration of ',...
                    'the RF pulse or increasing the TE and TR.'];

            else

                msg = ['The strength of the refocussing gradient pulse of the prepulse (',...
                    num2str(GZ_prepulse_ref_magn),' T/m)'...
                'is larger than the maximum gradient strength allowed (',...
                    num2str(G_max),' T/m) for the current configuration. Consider ',...
                    'increasing the TE and TR, decreasing the duration of ',...
                    'the RF pulse, increasing the maximum gradient strength allowed or ',...
                    'a combination of the above.'];

            end

            error(msg)

        end

        [GZ_prepulse_ref,GZ_prepulse_ref_plateau,GZ_prepulse_ref_area] = ...
            pulseSequenceGenerator.addGRtrapezoid(0,...
            (-1*GZ_prepulse_ref_magn),Gz_prepulse_ref_plateau_tmstps*dt,0,...
            'GZ_prepulse_ref',dt);

        block5      = zeros(8,N_block5);
        block5(6,:) = [GZ_prepulse,GZ_prepulse_ref];
        block5(1,GZ_prepulse_plateau(1,2):GZ_prepulse_plateau(1,3)) = RF2(1,:);
        block5(2,GZ_prepulse_plateau(1,2):GZ_prepulse_plateau(1,3)) = 90*pi/180; % phase 90 (in rads)
        
        % Update the timings for the other pulse sequence
        pulse_sequence  = [block5,pulse_sequence];

        isInKspace      = [zeros(1,size(block5,2)),isInKspace];

        kspace_times    = size(block5,2) + kspace_times;
        TR_times        = size(block5,2) + TR_times;

        soft_crushers   = [zeros(1,size(block5,2)),soft_crushers];
        
    elseif strcmp(bSSFPpreparation,'ramp')
        
        [pulse_sequence_RAMP,isInKspace_RAMP,~,...
            soft_crushers_RAMP,N_pulse_RAMP,~] = ...
            pulseSequenceGenerator.bSSFP_RampGenerator(RFmatrix,...
            acquiredkspace,FOV,receiver_BW,G_max,TR,TE,dt,gamma,...
            experiment_id,pulseq_id,conn_localdb,partialFourierStruct,...
            lengthRamp, GX, GX_pre, GX_ref, TEmin);
        
        % Update the timings for the other pulse sequence
        pulse_sequence  = [pulse_sequence_RAMP,pulse_sequence];

        isInKspace      = [isInKspace_RAMP,isInKspace];

        kspace_times    = N_pulse_RAMP + kspace_times;
        TR_times        = N_pulse_RAMP + TR_times;

        soft_crushers   = [soft_crushers_RAMP,soft_crushers];
        
    end

    
end

%% Add an a/2 at the end of bSSFP readout 
% This block will bring the magnetization back to the z axis. It consists
% of an extra TR
if addTailAby2
    
    % 1st block - Slice selection
    block1_tailAby2 = block1;

    % Since the phase of the first RF pulse of the readout is always -90,
    % the length_bSSFP will define whether the phase of the first RF pulse is
    % 90 or -90
    if mod(acquiredLines,2)==0  % even number 
        block1_tailAby2(2,GZ_plateau(1,2):GZ_plateau(1,3)) = -90*pi/180;
    elseif mod(acquiredLines,2)==1  % odd number
        block1_tailAby2(2,GZ_plateau(1,2):GZ_plateau(1,3)) = 90*pi/180;
    end

    % 2nd block - GRZ prephasor.
    block2_tailAby2 = zeros(8,size(GZ_ref,2));
    block2_tailAby2(6,1:size(GZ_ref,2)) = GZ_ref;

    % 3rd block - GRZ prephasor.
    block3_tailAby2      = block2_tailAby2;
    block3_tailAby2(6,:) = fliplr(block3_tailAby2(6,:));


    % 4th block - a2 RF and GRZ
    % Design RF
    [RFa2,~,~]    = pulseSequenceGenerator.addRFsinc(rf_dur,cycles,...
        angle/2,dt,gamma);

    block4_tailAby2 = zeros(8,max(size(RFa2,2),size(GZ,2)));
    block4_tailAby2(6,1:size(GZ,2)) = GZ;
    block4_tailAby2(1,GZ_plateau(1,2):GZ_plateau(1,3)) = RFa2(1,:);

    if mod(acquiredLines,2)==0  % even number 
        block4_tailAby2(2,GZ_plateau(1,2):GZ_plateau(1,3)) = 90*pi/180;
    elseif mod(acquiredLines,2)==1  % odd number
        block4_tailAby2(2,GZ_plateau(1,2):GZ_plateau(1,3)) = -90*pi/180;
    end
    
    %% Put all them together
    block_tailAby2 = zeros(8,N_TR);
    block_tailAby2(:,1:size(block1_tailAby2,2)) = block1_tailAby2;
    block_tailAby2(:,size(block1_tailAby2,2)+1:...
        size(block1_tailAby2,2)+size(block2_tailAby2,2)) = block2_tailAby2;
    block_tailAby2(:,floor(TE/dt)-size(block3_tailAby2,2)+1:floor(TE/dt)) = block3_tailAby2;
    block_tailAby2(:,floor(TE/dt)+1:floor(TE/dt)+size(block4_tailAby2,2)) = block4_tailAby2;

    soft_crushers_tailAby2   = zeros(1,size(block_tailAby2,2));
    isInKspace_tailAby2      = zeros(1,size(block_tailAby2,2));
    
    pulse_sequence  = [pulse_sequence,block_tailAby2];
    isInKspace      = [isInKspace,isInKspace_tailAby2];
    soft_crushers   = [soft_crushers,soft_crushers_tailAby2];
    
end

pulse_sequence(7,:) = 1:size(pulse_sequence,2);
pulse_sequence(8,:) = 1:size(pulse_sequence,2);
N_pulse             = size(pulse_sequence,2);

%% Calculate extra parameters
encodedkSpace.lines         = acquiredLines;
encodedkSpace.columns       = round(GX_plateau_duration/dt);

PulseqDuration = N_pulse*dt; %TR*acquiredLines;
disp(['The total duration of the pulse sequence is ',num2str(PulseqDuration),'sec'])
if mod(acquiredkspace(1,2),2)
    % If the number of kspace lines is an odd number, the center of
    % acquisition occurs at the center of the central kspace line.
    centerLine = round(acquiredkspace(1,2)/2);
    bSSFPCenterOfAcquisition = round((kspace_times(centerLine,2) + ...
        kspace_times(centerLine,1))/2);
else
    % If the number of kspace lines is an even number, the center of
    % acquisition occurs at the center between the centers of the lines 
    % that lie adjacent to the centre of kspace 
    centerLine1 = floor(acquiredkspace(1,2)/2);
    centerLine2 = ceil(acquiredkspace(1,2)/2);
    PulseqCenterOfAcquisition1 = round((kspace_times(centerLine1,2) + ...
        kspace_times(centerLine1,1))/2);
    PulseqCenterOfAcquisition2 = round((kspace_times(centerLine2,2) + ...
        kspace_times(centerLine2,1))/2);
    bSSFPCenterOfAcquisition = round((PulseqCenterOfAcquisition1 + ...
        PulseqCenterOfAcquisition2)/2);
end

%% FINAL CHECK
if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
        max(pulse_sequence(6,:))>=G_max
    exception = MException('GrMagn:GreaterThanMax',...
        'The gradients magnitude exceeds the maximum gradient strength');
    throw(exception)
end 