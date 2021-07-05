function [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,...
    N_TR,kspace_times,BW,TR_times,TRmin,PulseqDuration,encodedkSpace] = ...
    GREgenerator_OOP(RFmatrix,acquiredkspace,FOV,receiver_BW,G_max,TR,TE,dt,...
    gamma,conn_localdb,experiment_id,pulseq_id,partialFourierStruct,...
    deltaTime)
% This function generates GRE pulse sequences with in-phase and 
% out-of-phase imaging

N_TR            = floor(TR/dt);   % Number of time-steps during the TR

G_TR            = zeros(3,N_TR);
B1_TR           = zeros(2,N_TR);

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

% STEP 4 - GR_X pre.
% We assume that the pre-phasing pulse is a hard pulse with no ramps. The
% plateau duration is equal to the half duration of the "readout" 
% pulse and its magnitude is calculated so as the product to be equal to 
% the (GX_area(1,2)/2 + GX_area(1,1))

GX_pre_magn = 0.95*G_max;
if strcmp(partialFourierStruct.type,'readConjugate')
    GxFirstPart = (readConjugateFactor - 0.5) * (GX_area(1,2)/readConjugateFactor);
else
    GxFirstPart = GX_area(1,2)/2;
end
GX_pre_plateau_tmstps = floor((GxFirstPart + GX_area(1,1))/(GX_pre_magn*dt));
GX_pre_magn = (GxFirstPart + GX_area(1,1))/(GX_pre_plateau_tmstps*dt);
GX_pre_magn = GX_pre_magn*(-1);

[GX_pre,GX_pre_plateau,GX_pre_area] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GX_pre_magn,...
    GX_pre_plateau_tmstps*dt,0,'GX_pre',dt);

%% Calculate the min TE

block2_minDuration = (max(size(GZ_ref,2),size(GX_pre,2)))*dt;
TEmin = (size(block1,2)/2 + (readConjugateFactor - 0.5)*size(block3,2))*dt + block2_minDuration;

% Max. area of GY_PE
GY_ramp         = 0;
GY_max          = 0.95*G_max;
ky(1,1)         = ((acquiredkspace(1,2)-1)/2)*Dky;
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

block2_Duration = TE - (size(block1,2)/2 + (readConjugateFactor - 0.5)*size(block3,2))*dt;

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


%% 4th block - It holds the OOP readout gradients
% The prephasor of the second Greadout should have an area that 1) cancels 
% the area of the first readout after the formation of echo and 2) work as
% a prephasor for the second greadout
GxOOPFirstPart  = GxFirstPart;
GxSecondPart    = 0.5 * (GX_area(1,2)/readConjugateFactor);
GXOOP_pre_magn = 0.95*G_max;

GXOOP_pre_plateau_tmstps = floor((GxSecondPart+GX_area(1,3))+...
    (GxOOPFirstPart + GX_area(1,1))/(GXOOP_pre_magn*dt));
GXOOP_pre_magn = ((GxSecondPart+GX_area(1,3))+...
    (GxOOPFirstPart + GX_area(1,1)))/(GXOOP_pre_plateau_tmstps*dt);
GXOOP_pre_magn = GXOOP_pre_magn*(-1);

[GXOOP_pre,~,~] = ...
    pulseSequenceGenerator.addGRtrapezoid(0,GXOOP_pre_magn,...
    GXOOP_pre_plateau_tmstps*dt,0,'GX_pre',dt);

block4length    = round(deltaTime/dt);
block4minLength = size([GXOOP_pre,GX],2);

if block4length<block4minLength
    
    msg = ['The selected value of the time difference (delta) between the ',...
            'in-phase and out-of-phase echoes is too short. Please ',...
            'consider increasing it or decreasing the duration of the ',...
            'readout gradients.'];
        
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',...
        experiment_id,pulseq_id);
    
end

block4                              = zeros(8,block4length);
block4(4,end-block4minLength+1:end) = [GXOOP_pre,GX];


%% BUILD TR
TRmin = size(block1,2)+size(block2,2)+size(block3,2)+size(block4,2);

if N_TR<TRmin
    
    msg = ['The selected TR is lower that the minimum available TR (',...
            num2str(TRmin*dt),' sec) for the current configuration. Consider ',...
            'increasing the TR.'];
        
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
    
end

pulse_sequence_TR   = zeros(8,N_TR);
blocksLength        = size([block1,block2,block3,block4],2);

pulse_sequence_TR(:,1:blocksLength) = [block1,block2,block3,block4];

% pulse_sequence_TR(:,1:size(block1,2)) = block1;
% pulse_sequence_TR(:,size(block1,2)+1:...
%     size(block1,2)+size(block2,2)) = block2;
% pulse_sequence_TR(:,size(block1,2)+...
%     size(block2,2)+1:size(block1,2)+size(block2,2)+size(block3,2)) = block3;

soft_crushers_TR = zeros(1,N_TR);
soft_crushers_TR(1,N_TR) = 1;

pulse_sequence  = repmat(pulse_sequence_TR,1,acquiredLines);
soft_crushers   = repmat(soft_crushers_TR,1,acquiredLines);

isInKspace = zeros(1,size(pulse_sequence,2));

pulse_sequence(5,N_TR+1:end) = ...
    zeros(1,size(N_TR+1:size(pulse_sequence,2),2));

pulse_sequence(7,:) = 1:size(pulse_sequence,2);
pulse_sequence(8,:) = 1:size(pulse_sequence,2);

% kspace_times is a 3D matrix with two layers. The first layer holds the
% data for the fill on the in-phase kspace. The second layer holds the data
% for the fill of the out-of-phase kspace
kspace_times    = zeros(acquiredLines,2,2);

TR_times        = zeros(acquiredLines,2);

temp1 = 0;
for i=1:acquiredLines
    kspace_times(i,1,1) = ((i-1)*N_TR) + size(block1,2) + size(block2,2) + 1;
    kspace_times(i,2,1) = kspace_times(i,1,1) + (GX_plateau(1,3)-GX_plateau(1,2));
    
    kspace_times(i,1,2) = kspace_times(i,1,1) + block4length;
    kspace_times(i,2,2) = kspace_times(i,2,1) + block4length;
    
    TR_times(i,1) = (i-1)*N_TR+1;
    TR_times(i,2) = i*N_TR;
    
    isInKspace(1,kspace_times(i,1,1):kspace_times(i,2,1)) = temp1 + [1:GX_plateau(1,3)];
    temp1 = temp1 + GX_plateau(1,3);
    isInKspace(1,kspace_times(i,1,2):kspace_times(i,2,2)) = temp1 + [1:GX_plateau(1,3)];
    temp1 = temp1 + GX_plateau(1,3);
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
        GY_magn = GY_area(1,n)/block2_Duration;
    end
    
    [GY,GY_plateau,GY_area_2] = ...
        pulseSequenceGenerator.addGRtrapezoid(0,...
        GY_magn,block2_Duration,0,['GY',num2str(n)],dt);
                                           
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
    
end

%% Calculate extra parameters
encodedkSpace.lines         = acquiredLines;
encodedkSpace.columns       = round(GX_plateau_duration/dt);
encodedkSpace.contrasts     = 2;  % since we acquire two echoes per TR, one 
                                % in-phase and one out-of-phase

PulseqDuration          = TR*acquiredLines;
disp(['The total duration of the pulse sequence is ',num2str(PulseqDuration),'sec'])

%% FINAL CHECK
if max(pulse_sequence(4,:))>=G_max || max(pulse_sequence(5,:))>=G_max ||...
        max(pulse_sequence(6,:))>=G_max
    exception = MException('GrMagn:GreaterThanMax',...
        'The gradients magnitude exceeds the maximum gradient strength');
    throw(exception)
end