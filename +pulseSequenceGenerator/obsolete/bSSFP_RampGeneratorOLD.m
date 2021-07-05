function [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,...
    N_TR] = bSSFP_RampGenerator(RFmatrix,kspace,FOV,receiver_BW,G_max,...
    TR,TE,dt,gamma,lengthRamp)
% This function generates the ramp of bSSFP pulse sequences. 

N_TR            = floor(TR/dt);   % Number of time-steps during the TR

G_TR            = zeros(3,N_TR);
B1_TR           = zeros(2,N_TR);

rf_dur          = RFmatrix(1,1);
cycles          = RFmatrix(2,1);
angle           = RFmatrix(3,1);
slice_thickness = RFmatrix(4,1);

Dkx             = 1/FOV(1,1);   %k-space Dkx   ALLAXE TO
Dky             = 1/FOV(1,2);   %k-space Dky   ALLAXE TO

N_pulse         = kspace(1,2)*N_TR;
G_pulse         = zeros(3,N_pulse);
B1_pulse        = zeros(2,N_pulse);

flag            = 0;

pulse_sequence = [];

for iTR = 1:lengthRamp
    
    RFangle4ramp = (angle/(lengthRamp+1))*iTR;
    
    % 1st Block - RF and Gz (slice selection)

    % STEP 1 - Design RF
    [RF,BW,rf_timesteps]    = pulseSequenceGenerator.addRFsinc(rf_dur,cycles,...
        RFangle4ramp,dt,gamma);

    % STEP 2 - Design slice selection gradient (GZ)
    Gz_magn                 = BW/(gamma*slice_thickness);
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
    
    % Since the phase of the first RF pulse of the readout is always -90,
    % the lengthRamp will define whether the phase of the first RF pulse of
    % the ramp is 90 or -90
    if mod(lengthRamp,2)==0  % even number 
        sign = cosd(iTR*180);
        block1(2,GZ_plateau(1,2):GZ_plateau(1,3)) = sign*90*pi/180;
    elseif mod(lengthRamp,2)==1  % odd number
        sign = cosd(iTR*180);
        block1(2,GZ_plateau(1,2):GZ_plateau(1,3)) = -sign*90*pi/180;
    end

    % 3rd Block - GR_X readout

    % STEP 5 - GR_X
    GX_magn                 = receiver_BW/(gamma*FOV(1,1));
    GX_plateau_duration     = kspace(1,1)/receiver_BW;  
    [GX,GX_plateau,GX_area] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,GX_magn,...
        GX_plateau_duration,0,'GX',dt);

    block3 = zeros(8,size(GX,2));
    block3(4,1:size(GX,2)) = GX;

    % 2nd Block - GR_Z ref. + GR_X pref.

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

    % Calculate the min TE

    block2_minDuration = (max(size(GZ_ref,2),size(GX_pre,2)))*dt;
    TEmin = (size(block1,2)/2 + size(block3,2)/2)*dt + block2_minDuration;

    % Max. area of GY_PE
    GY_ramp         = 0;
    GY_max          = 0.95*G_max;
    ky(1,1)         = ((kspace(1,2)-1)/2)*Dky;
    GY_area(1,1) 	= ky(1,1)/gamma; 

    while GY_area(1,1)/block2_minDuration > GY_max

        block2_minDuration  = block2_minDuration + dt;
        TEmin               = TEmin + dt;
        flag                = 1;

    end

    if TE < TEmin

        if flag == 1

            msg = ['The selected TE is smaller than the minimum available TE (',...
                num2str(TEmin),'sec) for the current configuration. Consider',...
                ' incresing the maximum gradient strength.'];        

        else

            msg = ['The selected TE is smaller than the minimum available TE (',...
                num2str(TEmin),'sec) for the current configuration. Consider',...
                ' decreasing the duration of the RF pulse or the duration of',...
                ' the readout.'];

        end

        error(msg)

    end

    % 2nd Block - PE

    block2_Duration = TE - (size(block1,2)/2 + size(block3,2)/2)*dt;

    % STEP 6 - PE
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

    block2 = zeros(8,floor(block2_Duration/dt));
    block2(6,1:size(GZ_ref,2)) = GZ_ref;
    block2(4,end-size(GX_pre,2)+1:end) = GX_pre;
    block2(5,:) = GY;

    % 4th block - GRY refoc. + PE refoc. + GRZ refoc.
    block4      = block2;
    block4(5,:) = -block4(5,:);
    block4(6,:) = fliplr(block4(6,:));

    % BUILD TR
    TRmin = size(block1,2)+2*size(block2,2)+size(block3,2);

    if N_TR<TRmin

        msg = ['The selected TR is lower that the minimum available TR (',...
                num2str(TRmin),' sec) for the current configuration. Consider ',...
                'increasing the TR.'];

        error(msg)

    end

    pulse_sequence_TR = zeros(8,N_TR);
    pulse_sequence_TR(:,1:size(block1,2)) = block1;
    pulse_sequence_TR(:,size(block1,2)+1:...
        size(block1,2)+size(block2,2)) = block2;
    pulse_sequence_TR(:,size(block1,2)+...
        size(block2,2)+1:size(block1,2)+size(block2,2)+size(block3,2)) = block3;
    pulse_sequence_TR(:,size(block1,2)+...
        size(block2,2)+size(block3,2)+1:size(block1,2)+size(block2,2)+...
        size(block3,2)+size(block4,2)) = block4;
    
    pulse_sequence = [pulse_sequence,pulse_sequence_TR];
end

soft_crushers   = zeros(1,size(pulse_sequence,2));
isInKspace      = zeros(1,size(pulse_sequence,2));

pulse_sequence(7,:) = 1:size(pulse_sequence,2);
pulse_sequence(8,:) = 1:size(pulse_sequence,2);
N_pulse             = size(pulse_sequence,2);

%% FINAL CHECK
if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
        max(pulse_sequence(6,:))>=G_max
    exception = MException('GrMagn:GreaterThanMax',...
        'The gradients magnitude exceeds the maximum gradient strength');
    throw(exception)
end