function [pulse_sequence,isInKspace,soft_crushers,N_pulse,info,dt,gamma,...
    tirl,printForImage,elementsize_new,slabThickness,dwPulseSequence] = ...
    createPulseSequence(pulseSeqFamilyName,struct_pulseq,conn_localdb,...
    performance_gridZ_sliceThickness,dwell_time,experiment_id,pulseq_id,...
    noOfSlices,advancedNotifications)
% printForImage holds the information written on top of the viewer. This
% information is related to the pulse sequence

dwPulseSequence = [];

% Check if the selected value of the receiver BW is included in the list of
% available values that will not create artifacts due to rounding errors.
% The rounding errors arise from the fact that the duration of the readout
% gradient is not an integer multiple of dt. The check is performed for
% dt=10^-6. The duration of the readout gradient is defined by the 
% kspace(1,1)/receiver_BW
receiverBW  = str2num(struct_pulseq.BW);
kspace(1,1)	= str2num(struct_pulseq.matrix_x);  % defines the size of k-space along x direction
kspace(1,2)	= str2num(struct_pulseq.matrix_y);  % defines the size of k-space along y direction
dt_test         = 10^-6;
rBW_allowList   = [];
for rBW_test=50000:250000
    if rem(kspace(1,1)/rBW_test,dt_test)==0
        rBW_allowList = [rBW_allowList,rBW_test];
    end
end
if ~ismember(receiverBW,rBW_allowList)
    msg = ['The selected receiver BW is not allowed since it ',...
        'may cause artifacts due to rounding errors.'];
    eduTool.frontend.errorAndDBentry(conn_localdb,msg,'cancelled-error',...
        experiment_id,pulseq_id);
end
% If the user has selected a fixed dt, then dt is equal to 10^-6. If the
% user has selected a dynamic dt, dt will be defined by the selected
% receiver BW so as the duration of the readout gradient to be an integer 
% multiple of dt (in order to avoid artifacts due to rounding errors). The 
% duration of the readout gradient is defined by the kspace(1,1)/receiver_BW
if strcmp(dwell_time,'fixed')
    dt = 10^-6;
elseif strcmp(dwell_time,'dynamic')
    dt_allowList = [];
    for i=1:0.1:5
        dt_new = i*10^-6;
        if rem(kspace(1,1)/receiverBW,dt_new)==0
            dt_allowList = [dt_allowList,dt_new];
        end
    end
%     if receiverBW<=100000 
%         dt = max(dt_allowList);
%     elseif receiverBW>100000 && receiverBW<=200000
%         dt_options = dt_allowList(find(dt_allowList<=2.5*10^-6));
%         dt = max(dt_options);
%     elseif receiverBW>200000 && receiverBW<=250000
%         dt_options = dt_allowList(find(dt_allowList<=2*10^-6));
%         dt = max(dt_options);
%     else
    if receiverBW == 50000
        dt = 5*10^-6;
    elseif receiverBW == 100000
        dt = 5*10^-6;
    elseif receiverBW == 200000
        dt = 2.5*10^-6;
    elseif receiverBW == 250000
        dt = 2*10^-6;
    else
        dt = 10^-6;
    end
end
disp(['dt = ',num2str(dt),'sec'])

FOV = [str2num(struct_pulseq.fov_x),str2num(struct_pulseq.fov_y)];

acquiredkspace  = kspace;
acquiredFOV     = FOV;

if isfield(struct_pulseq,'foldoversuppr')
    if strcmp(struct_pulseq.foldoversuppr,'yes')
        acquiredkspace(1,2) = 2*kspace(1,2);
        acquiredFOV(1,2)    = 2*FOV(1,2);
    end
end

if isfield(struct_pulseq,'partialFourier')
    partialFourierStruct.type = struct_pulseq.partialFourier;
    if ~strcmp(partialFourierStruct.type,'no')
        partialFourierStruct.factor = ...
            str2num(struct_pulseq.partialFourierFactor);
        % If read-conjugate active, adjust the partialFourierFactor so as
        % to result in an integer number of Greadout timesteps
        if strcmp(partialFourierStruct.type,'readConjugate') && ...
                ~isequal(partialFourierStruct.factor,1)
            GreadoutDuration = partialFourierStruct.factor*...
                kspace(1,1)/receiverBW;
            if rem(GreadoutDuration,dt)~=0
                GreadoutDuration_tmstps = GreadoutDuration/dt;
                partialFourierStruct.factor = round(GreadoutDuration_tmstps)/...
                    ((kspace(1,1)/receiverBW)/dt);
            end
        end
    else
        partialFourierStruct.factor = 1;  % in this case, no half-fourier is applied
    end
    % Create a text string for printing this on the image viewers
    if strcmp(partialFourierStruct.type,'readConjugate')
        partialFourierText = [' (',num2str(round(partialFourierStruct.factor*kspace(1,1))),...
            'x',num2str(kspace(1,2)),')'];
    elseif strcmp(partialFourierStruct.type,'phaseConjugate')
        partialFourierText = [' (',num2str(kspace(1,1)),...
            'x',num2str(round(partialFourierStruct.factor*kspace(1,2))),')'];
    else
        partialFourierText = '';
    end
else
    partialFourierStruct.type   = 'no';
    partialFourierStruct.factor = 1;  % in this case, no half-fourier is applied
    partialFourierText = '';
end

% Use parallel imaging, only when no half-fourier is applied
if isfield(struct_pulseq,'parallelImaging') ...
        && strcmp(struct_pulseq.parallelImaging,'sense')
    Rfactor             = str2num(struct_pulseq.rfactor);
    if rem(kspace(1,2),Rfactor)~=0
        msg = ['The selected size of the image along the PE direction ',...
            '(Matrix - Phase) is not an integer multiple of the parallel ',...
            'imaging acceleration factor (R-factor). Please modify either the ',...
            'size of the image along the PE direction or the acceleration ',...
            'factor (R-factor).'];
        eduTool.frontend.errorAndDBentry(conn_localdb,msg,'cancelled-error',...
            experiment_id,pulseq_id);
    else
        acquiredkspace(1,2) = kspace(1,2)/Rfactor;
        acquiredFOV(1,2)    = FOV(1,2)/Rfactor;
    end
end

if strcmp(pulseSeqFamilyName,'GRE')
    
    TR              = str2num(struct_pulseq.tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);

    fast = 0;  % CHANGE IT
    
    gamma           = 42.56*10^6;                           % Hz/T
    f0              = 64000000;                             % in Hz

    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = sliceThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m
    %dt              = 4*10^-6;%1/receiverBW;%10^-6;                                % sec
    
    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    slabThickness   = sliceThickness;
    
    [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
        kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
        pulseSequenceGenerator.GREgenerator_PerformanceOpt(RFmatrix,acquiredkspace,acquiredFOV,...
        receiverBW,G_max,TR,TE,dt,gamma,conn_localdb,experiment_id,pulseq_id,...
        partialFourierStruct);
    info = pulseSequenceGenerator.createInfoStruct_GRE(kspace,FOV,dt,...
        receiverBW,TR,TE,kspace_times,f0,sliceThickness,acquiredkspace,...
        acquiredFOV,encodedkSpace);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_GRE(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,acquiredkspace,...
            acquiredFOV,encodedkSpace);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2)),' ',partialFourierText];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
    else
        printForImage.fourthLine    = '';
    end
elseif strcmp(pulseSeqFamilyName,'GRE-conc')    
    % GRE with concatenations        
    
    TR              = str2num(struct_pulseq.tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    sliceGap        = str2num(struct_pulseq.slice_gap);
    concatenations  = str2num(struct_pulseq.concatenations);
    slices          = str2num(struct_pulseq.slices);
    
    kspace(1,3)     = str2num(struct_pulseq.slices);
    acquiredkspace(1,3) = kspace(1,3);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    f0              = 64000000;                             % in Hz

    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = sliceThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m
    %dt              = 4*10^-6;%1/receiverBW;%10^-6;                                % sec
    
    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    % Important for a 3D pulse sequence not to be forgotten
    slabThickness   = kspace(1,3)*sliceThickness+(kspace(1,3)-1)*sliceGap;
    
    [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
        kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace,RFfreqAmpl] = ...
        pulseSequenceGenerator.GREgenerator_concatenations(RFmatrix,...
        acquiredkspace,acquiredFOV,receiverBW,G_max,TR,TE,dt,gamma,...
        conn_localdb,experiment_id,pulseq_id,partialFourierStruct,...
        slices,concatenations,sliceGap);
    
%     [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
%         kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
%         pulseSequenceGenerator.GREgenerator(RFmatrix,acquiredkspace,acquiredFOV,...
%         receiverBW,G_max,TR,TE,dt,gamma,conn_localdb,experiment_id,pulseq_id,...
%         partialFourierStruct);
    info = pulseSequenceGenerator.createInfoStruct_GRE_conc(kspace,FOV,dt,...
        receiverBW,TR,TE,kspace_times,f0,sliceThickness,acquiredkspace,...
        acquiredFOV,encodedkSpace,RFfreqAmpl);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_GRE_conc(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,acquiredkspace,...
            acquiredFOV,encodedkSpace,RFfreqAmpl);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2)),' ',partialFourierText];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
    else
        printForImage.fourthLine    = '';
    end
elseif strcmp(pulseSeqFamilyName,'GRE-OOP')
    % GRE out-of-phase. For more info, check pages 860-862 from Handbook of
    % MRI pulse Sequences.
    
    TR              = str2num(struct_pulseq.tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    if isfield(struct_pulseq,'OOP_deltaTime')
        deltaTime = str2num(struct_pulseq.OOP_deltaTime);
    else
        msg = ['The value of the time difference (delta) between the ',...
            'in-phase and out-of-phase echoes is missing. Please contact ',...
            'the administrator.'];
        eduTool.frontend.errorAndDBentry(conn_localdb,msg,'cancelled-error',...
            experiment_id,pulseq_id);
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    f0              = 64000000;                             % in Hz

    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = sliceThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m
    %dt              = 4*10^-6;%1/receiverBW;%10^-6;                                % sec
    
    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    slabThickness   = sliceThickness;
    
    [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
        kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
        pulseSequenceGenerator.GREgenerator_OOP(RFmatrix,acquiredkspace,acquiredFOV,...
        receiverBW,G_max,TR,TE,dt,gamma,conn_localdb,experiment_id,pulseq_id,...
        partialFourierStruct,deltaTime);
    info = pulseSequenceGenerator.createInfoStruct_GRE(kspace,FOV,dt,...
        receiverBW,TR,TE,kspace_times,f0,sliceThickness,acquiredkspace,...
        acquiredFOV,encodedkSpace);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_GRE(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,acquiredkspace,...
            acquiredFOV,encodedkSpace);
    end
    
    % We are expecting 4 images: 1 simulated image from 1st echo, 1
    % simulated image from 2nd echo, the water image and the fat image
    info.reconstruction.imagesPerSlice = 4;
    
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2)),' ',partialFourierText];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
    else
        printForImage.fourthLine    = '';
    end
elseif strcmp(pulseSeqFamilyName,'GRE-3D')
    
    kspace(1,3) = str2num(struct_pulseq.slices);
    acquiredkspace(1,3) = kspace(1,3);
    
    TR              = str2num(struct_pulseq.tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    f0              = 64000000;                             % in Hz

    slabThickness   = kspace(1,3)*sliceThickness;
    acquiredFOV(1,3)= slabThickness;
        
    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = slabThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m
   
    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
        kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
        pulseSequenceGenerator.GREgenerator_3D(RFmatrix,...
        kspace,FOV,receiverBW,G_max,TR,TE,dt,gamma,conn_localdb,...
        experiment_id,pulseq_id);
    info = pulseSequenceGenerator.createInfoStruct_GRE_3D(kspace,FOV,dt,receiverBW,TR,TE,...
        kspace_times,f0,sliceThickness,acquiredkspace,acquiredFOV,...
        encodedkSpace);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_GRE_3D(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,acquiredkspace,...
            acquiredFOV,encodedkSpace);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];
    printForImage.fourthLine    = '';
    
elseif strcmp(pulseSeqFamilyName,'MP-RAGE')
    
    kspace(1,3) = str2num(struct_pulseq.slices);
    acquiredkspace(1,3) = kspace(1,3);
    
    TR              = str2num(struct_pulseq.tr);
    MPRAGETR        = str2num(struct_pulseq.mp_rage_tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    segments        = kspace(1,2);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    f0              = 64000000;                             % in Hz

    slabThickness   = kspace(1,3)*sliceThickness;
    acquiredFOV(1,3)= slabThickness;
        
    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = slabThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m
   
    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
        kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
        pulseSequenceGenerator.MPRAGEgenerator_GRE3D(RFmatrix,...
        kspace,FOV,receiverBW,G_max,TR,TE,dt,gamma,MPRAGETR,...
        struct_pulseq,segments,conn_localdb,experiment_id,pulseq_id);
    info = pulseSequenceGenerator.createInfoStruct_GRE_3D(kspace,FOV,dt,receiverBW,TR,TE,...
        kspace_times,f0,sliceThickness,acquiredkspace,acquiredFOV,...
        encodedkSpace);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_GRE_3D(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,acquiredkspace,...
            acquiredFOV,encodedkSpace);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];
    printForImage.fourthLine    = '';
    
elseif strcmp(pulseSeqFamilyName,'EPI')
    
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    %dt              = 10^-6;                                % sec
    f0              = 64000000;                             % in Hz

    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = sliceThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;
    
    slabThickness   = sliceThickness;
    
    % The halfFourierFactor should give you an even number of acquired
    % kspace lines
    acquiredLines   = round(partialFourierStruct.factor*acquiredkspace(1,2));
    if rem(acquiredLines,2)==1
        partialFourierStruct.factor = (acquiredLines+1)/acquiredkspace(1,2);
    end
    
    [pulse_sequence,isInKspace,soft_crushers,N_pulse,...
        kspace_times,BW_rf,PulseqDuration,ESP,encodedkSpace] = ...
        pulseSequenceGenerator.EPIgenerator(RFmatrix,acquiredkspace,...
        acquiredFOV,receiverBW,G_max,dt,gamma,conn_localdb,experiment_id,...
        pulseq_id,partialFourierStruct);
    
    info = pulseSequenceGenerator.createInfoStruct_EPI(kspace,FOV,dt,...
        receiverBW,kspace_times,f0,sliceThickness,ESP,acquiredkspace,...
        acquiredFOV,partialFourierStruct,encodedkSpace);
    
    tirl                = N_pulse*dt;
    TRmin               = N_pulse;
    pulse_sequence_TR   = pulse_sequence;
    %% ADD FAST ALGORITHM
    if fast == 1
        if advancedNotifications
            eduTool.frontend.errorAndDBentry(conn_localdb,...
                'The fast algorithm has been deactivated for this pulse sequence.',...
                'info',experiment_id,pulseq_id);
        end
%         [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
%             pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
%             info.pulseSequence.kspace,soft_crushers);    
%         info = pulseSequenceGenerator.createInfoStruct_EPI(kspace,FOV,dt,...
%             receiverBW,kspace_times_new,f0,sliceThickness,ESP);
%         N_pulse = size(pulse_sequence,2);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];
    printForImage.thirdLine     = ['ESP: ',num2str(ESP*1000),'ms'];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
    else
        printForImage.fourthLine    = '';
    end
elseif strcmp(pulseSeqFamilyName,'SE-EPI') 
    
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    f0              = 64000000;                             % in Hz

    RFmatrix1        = zeros(6,1);
    RFmatrix1(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix1(2,:)   = 2;                                    % RF sinc cycles
    if isfield(struct_pulseq,'RFflipAngle')
        RFmatrix1(3,:)   = str2num(...
            struct_pulseq.RFflipAngle);
    else
        RFmatrix1(3,:)   = 90;          % RF rotation angle
    end
    RFmatrix1(4,:)   = sliceThickness;                       % desired slice thickness
    RFmatrix1(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix1(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis
    
    RFmatrix2        = zeros(6,1);
    RFmatrix2(1,:)   = str2num(...
        struct_pulseq.RFduration);       % RF duration 
    RFmatrix2(2,:)   = 2;                % RF sinc cycles
    if isfield(struct_pulseq,'refoc_pulse_angle')
        RFmatrix2(3,:)   = str2num(...
            struct_pulseq.refoc_pulse_angle);
    else
        RFmatrix2(3,:)   = 180;          % RF rotation angle
    end
    RFmatrix2(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
    RFmatrix2(5,:)   = [1];              % starting timestep of RF
    RFmatrix2(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis
    
    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m

    TE              = str2num(struct_pulseq.te);
    
    slabThickness   = sliceThickness;
    
    % The halfFourierFactor should give you an even number of acquired
    % kspace lines
    acquiredLines   = round(partialFourierStruct.factor*acquiredkspace(1,2));
    if rem(acquiredLines,2)==1
        partialFourierStruct.factor = (acquiredLines+1)/acquiredkspace(1,2);
    end
    
    [pulse_sequence,isInKspace,soft_crushers,N_pulse,...
        kspace_times,BW_rf,PulseqDuration,ESP,encodedkSpace] = ...
        pulseSequenceGenerator.EPI_SEgenerator(RFmatrix1,RFmatrix2,...
        acquiredkspace,acquiredFOV,receiverBW,...
        G_max,dt,TE,gamma,conn_localdb,experiment_id,pulseq_id,...
        partialFourierStruct);
    
    info = pulseSequenceGenerator.createInfoStruct_EPI(kspace,FOV,dt,...
        receiverBW,kspace_times,f0,sliceThickness,ESP,acquiredkspace,...
        acquiredFOV,partialFourierStruct,encodedkSpace);
    
    tirl                = N_pulse*dt;
    TRmin               = N_pulse;
    pulse_sequence_TR   = pulse_sequence;
    %% ADD FAST ALGORITHM
    if fast == 1
        if advancedNotifications
            eduTool.frontend.errorAndDBentry(conn_localdb,...
                'The fast algorithm has been deactivated for this pulse sequence.',...
                'info',experiment_id,pulseq_id);
        end
%         [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
%             pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
%             info.pulseSequence.kspace,soft_crushers);    
%         info = pulseSequenceGenerator.createInfoStruct_EPI(kspace,FOV,dt,...
%             receiverBW,kspace_times_new,f0,sliceThickness,ESP);
%         N_pulse = size(pulse_sequence,2);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];
    printForImage.thirdLine     = ['ESP: ',num2str(ESP*1000),'ms'];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
    else
        printForImage.fourthLine    = '';
    end
    
elseif strcmp(pulseSeqFamilyName,'bSSFP') || strcmp(pulseSeqFamilyName,'IR-bSSFP')
    
    TR              = str2num(struct_pulseq.tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    gamma           = 42.56*10^6;               % Hz/T
    %dt              = 10^-6;                    % sec
    f0              = 64000000;                 % in Hz

    RFmatrix        = zeros(6,1);
    RFmatrix(1,:)   = str2num(struct_pulseq.RFduration);    % RF duration 
    RFmatrix(2,:)   = 2;                                    % RF sinc cycles
    RFmatrix(3,:)   = str2num(struct_pulseq.RFflipAngle);   % RF rotation angle
    RFmatrix(4,:)   = sliceThickness;                       % desired slice thickness
    RFmatrix(5,:)   = 1;                                    % starting timestep of RF
    RFmatrix(6,:)   = 1;                                    % 1 is for x-axis, 2 is for y-axis

    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m

    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    slabThickness   = sliceThickness;
    if ~strcmp(pulseSeqFamilyName,'MOLLI')
        [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
            kspace_times,BW_rf,TR_times,receiver_phase,TRmin,tirl,...
            bSSFPCenterOfAcquisition,encodedkSpace] = ...
            pulseSequenceGenerator.bSSFPgenerator(RFmatrix,...
            acquiredkspace,acquiredFOV,receiverBW,G_max,TR,TE,dt,gamma,...
            struct_pulseq.bSSFPpreparation,experiment_id,pulseq_id,conn_localdb,...
            partialFourierStruct,str2num(struct_pulseq.bSSFPrampLength),0);

        if (strcmp(pulseSeqFamilyName,'IR-bSSFP') && isfield(struct_pulseq,'ir_time')) || ...
                (isfield(struct_pulseq,'fatsat') && strcmp(struct_pulseq.fatsat,'lipidir'))
            if 1
                
                % In IR-bSSFP the inversion time counts from the end of the IR 
                % pulse till the middle of the central kspace line of the bSSFP 
                % acquisition.
                TIendPoint          = bSSFPCenterOfAcquisition;

                structIR.TI         = str2num(struct_pulseq.ir_time);
                if isfield(struct_pulseq,'fatsat') && strcmp(struct_pulseq.fatsat,'lipidir')
                    structIR.type       = 'sinc';
                    structIR.IRduration = 0.003;
                else
                    structIR.type       = struct_pulseq.ir_type;
                    structIR.IRduration = str2num(struct_pulseq.ir_duration);
                end
                structIR.angle      = 180;
                structIR.phase      = 0; % in degrees
                structIR.IRcycles   = 2; % only for sinc RF pulse

                structExper.dt      = dt;
                structExper.gamma   = gamma;

                compressIRmodule    = 1;

                % blockIR_actual_timesteps holds the "true" number of timesteps of the
                % blockIR. blockIR_timesteps holds the "effective" number of timesteps due
                % to the application of the fast algorithm

                [blockIR,blockIR_actual_timesteps,blockIR_timesteps] = ...
                    pulseSequenceGenerator.generateIRblock(structIR,...
                    structExper,TIendPoint,'bSSFP',compressIRmodule,conn_localdb,...
                    experiment_id,pulseq_id);
                
            else
                
                % In IR-bSSFP the inversion time counts from the end of the IR 
                % pulse till the middle of the central kspace line of the bSSFP 
                % acquisition.
                TI = str2num(struct_pulseq.ir_time); %0.69*0.312 - ln(2) * T1fat
                if isfield(struct_pulseq,'fatsat') && strcmp(struct_pulseq.fatsat,'lipidir')
                    IRtype      = 'sinc';
                    IRduration  = 0.003;
                elseif (strcmp(pulseSeqFamilyName,'IR-bSSFP') && isfield(struct_pulseq,'ir_time'))
                    IRtype      = struct_pulseq.ir_type;
                    IRduration  = str2num(struct_pulseq.ir_duration);
                end

                TIendPoint = bSSFPCenterOfAcquisition;

                [blockIR,blockIR_timesteps,blockIR_actual_timesteps] = ...
                    pulseSequenceGenerator.IRgenerator(dt,gamma,TI,...
                    IRduration,IRtype,TIendPoint,conn_localdb,...
                    experiment_id,pulseq_id,pulseSeqFamilyName);
                
            end

            pulse_sequence(7,:) = pulse_sequence(7,:) + blockIR_timesteps;
            pulse_sequence(8,:) = pulse_sequence(8,:) + blockIR_timesteps;
            pulse_sequence = [blockIR,pulse_sequence];
            isInKspace = [zeros(1,size(blockIR,2)),isInKspace];
            soft_crushers = [zeros(1,size(blockIR,2)),soft_crushers];

            kspace_times = kspace_times + blockIR_timesteps;

            N_pulse = size(pulse_sequence,2);

            tirl = tirl + blockIR_actual_timesteps * dt;

        end

        info = pulseSequenceGenerator.createInfoStruct_bSSFP(kspace,FOV,dt,receiverBW,TR,TE,...
            kspace_times,f0,sliceThickness,receiver_phase,RFmatrix,acquiredkspace,...
            acquiredFOV,encodedkSpace);
        
    else
        % Create first the bSSFP readout. It keeps the preparation (ramp or
        % a/2), the actual bSSFP and the tail a/2.
        [pulse_sequence_bSSFP,isInKspace_bSSFP,pulse_sequence_TR,...
            soft_crushers_bSSFP,N_pulse_bSSFP,N_TR_bSSFP,...
            kspace_times_bSSFP,BW_rf,TR_times_bSSFP,receiver_phase_bSSFP,...
            TRmin,tirl_bSSFP,bSSFPCenterOfAcquisition,encodedkSpace] = ...
            pulseSequenceGenerator.bSSFPgenerator(RFmatrix,...
            acquiredkspace,acquiredFOV,receiverBW,G_max,TR,TE,dt,gamma,...
            struct_pulseq.bSSFPpreparation,experiment_id,pulseq_id,conn_localdb,...
            partialFourierStruct,str2num(struct_pulseq.bSSFPrampLength),1);
        
        % Create the info structure
        info_bSSFP = pulseSequenceGenerator.createInfoStruct_bSSFP(kspace,FOV,dt,receiverBW,TR,TE,...
            kspace_times_bSSFP,f0,sliceThickness,receiver_phase_bSSFP,RFmatrix,acquiredkspace,...
            acquiredFOV,encodedkSpace);
        
        pulse_sequence_READOUT.pulse_sequence   = pulse_sequence_bSSFP;
        pulse_sequence_READOUT.isInKspace       = isInKspace_bSSFP;
        pulse_sequence_READOUT.soft_crushers    = soft_crushers_bSSFP;
        pulse_sequence_READOUT.N_pulse          = N_pulse_bSSFP;
        pulse_sequence_READOUT.kspace_times     = kspace_times_bSSFP;
        pulse_sequence_READOUT.TR_times         = TR_times_bSSFP;
        pulse_sequence_READOUT.receiver_phase   = receiver_phase_bSSFP;
        
        % Formulate the MOLLI pulse sequence
        cardiac_cycle_duration_ms   = 1000; % in msec
        cardiac_cycle_duration_s    = cardiac_cycle_duration_ms/1000; % in sec
        [MOLLI_scheme,pause_cc,TI_initial] = ...
            pulseSequenceGenerator.extractMOLLIspecs(struct_pulseq,cardiac_cycle_duration_ms);
        
        % Check if the TIs have not been given in an ascending order
        if TI_initial~=sort(TI_initial)
            msg = ['The TIs are not in an ascending order. Please modify the TIs of the selected MOLLI scheme.'];

            eduTool.frontend.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
        end

        TD_initial              = 0.5;  % choose a low TD_initial to fit the bSSFP readout 
                                        % within the cardiac cycle
        
        % TD is the time from the R-wave triggering till the start of the ramp
        TD                  = abs(repmat(TD_initial,1,size(TI_initial,2)));
        
        % Inversion pulse%         IRduration              = 0.01024; % for HypSec
%         IRduration              = 0.00256; % for TanTanh
%         IRduration              = 0.005; % for BIR4
%         structIR.type       = 'TanTanh';%'BIR4'; % %'HS1'; % 'TanTanh'
%         structIR.IRduration = IRduration;
        
        structIR.type       = struct_pulseq.ir_type_molli;
        structIR.IRduration = str2num(struct_pulseq.ir_duration_molli);
        structIR.angle      = 180;
        
        structExper.dt      = dt;
        structExper.gamma   = gamma;
        
        blockIR = pulseSequenceGenerator.generateIRpulse(structIR,structExper);
        
        % Create MOLLI pulse sequence
        if fast==1
            [pulse_sequence_MOLLI,times_fitting_single_point,tirl] = ...
                pulseSequenceGenerator.MOLLIgenerator_ImprovedPerformance(pulse_sequence_READOUT,...
                info_bSSFP,blockIR,cardiac_cycle_duration_s,TD,...
                TI_initial,dt,MOLLI_scheme,pause_cc,TE,TR,RFmatrix(1,:),acquiredkspace,...
                conn_localdb,experiment_id,pulseq_id);
        else
            [pulse_sequence_MOLLI,times_fitting_single_point,tirl] = ...
                pulseSequenceGenerator.MOLLIgenerator(pulse_sequence_READOUT,...
                info_bSSFP,blockIR,cardiac_cycle_duration_s,TD,...
                TI_initial,dt,MOLLI_scheme,pause_cc,TE,TR,RFmatrix(1,:),acquiredkspace,...
                conn_localdb,experiment_id,pulseq_id);
        end
        
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
        
    end
    tic
    if fast == 1 && ~strcmp(pulseSeqFamilyName,'MOLLI')
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_bSSFP(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,receiver_phase,...
            RFmatrix,acquiredkspace,acquiredFOV,encodedkSpace);
    end
    toc
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];    
    if strcmp(struct_pulseq.parallelImaging,'sense')
        if (strcmp(pulseSeqFamilyName,'IR-bSSFP') && isfield(struct_pulseq,'ir_time')) || ...
                (isfield(struct_pulseq,'fatsat') && strcmp(struct_pulseq.fatsat,'lipidir'))
            printForImage.fourthLine = ['TI: ',num2str(str2num(struct_pulseq.ir_time)*1000),'ms',' - SENSE: ',struct_pulseq.rfactor];
        else
            printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
        end
    else
        if (strcmp(pulseSeqFamilyName,'IR-bSSFP') && isfield(struct_pulseq,'ir_time')) || ...
                (isfield(struct_pulseq,'fatsat') && strcmp(struct_pulseq.fatsat,'lipidir'))
            printForImage.fourthLine = ['TI: ',num2str(str2num(struct_pulseq.ir_time)*1000),'ms'];
        else
            printForImage.fourthLine = '';
        end
    end    

elseif strcmp(pulseSeqFamilyName,'SE') || strcmp(pulseSeqFamilyName,'IR-SE')
    
    TR              = str2num(struct_pulseq.tr);
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    %dt              = 10^-6;                                % sec
    f0              = 64000000;                             % in Hz

    RFmatrix1        = zeros(6,1);
    RFmatrix1(1,:)   = str2num(...
        struct_pulseq.RFduration);       % RF duration
    RFmatrix1(2,:)   = 2;                % RF sinc cycles
    RFmatrix1(3,:)   = 90;               % RF rotation angle
    RFmatrix1(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
    RFmatrix1(5,:)   = [1];              % starting timestep of RF
    RFmatrix1(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis

    RFmatrix2        = zeros(6,1);
    RFmatrix2(1,:)   = str2num(...
        struct_pulseq.RFduration);       % RF duration 
    RFmatrix2(2,:)   = 2;                % RF sinc cycles
    if isfield(struct_pulseq,'refoc_pulse_angle')
        RFmatrix2(3,:)   = str2num(...
            struct_pulseq.refoc_pulse_angle);
    else
        RFmatrix2(3,:)   = 180;          % RF rotation angle
    end
    RFmatrix2(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
    RFmatrix2(5,:)   = [1];              % starting timestep of RF
    RFmatrix2(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis
    
    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m

    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    slabThickness   = sliceThickness;
    
    if strcmp(pulseSeqFamilyName,'IR-SE')
        TI          = str2num(struct_pulseq.ir_time);
        IRtype      = struct_pulseq.ir_type;
        IRduration  = str2num(struct_pulseq.ir_duration);
        
        [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
            kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
            pulseSequenceGenerator.SEgeneratorWithIRoption(RFmatrix1,RFmatrix2,...
            acquiredkspace,acquiredFOV,receiverBW,G_max,TR,TE,dt,gamma,...
            conn_localdb,experiment_id,pulseq_id,partialFourierStruct,1,TI,IRtype,...
            IRduration,struct_pulseq);
        
    elseif strcmp(pulseSeqFamilyName,'SE')
        
        [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
            kspace_times,BW_rf,TR_times,TRmin,tirl,encodedkSpace] = ...
            pulseSequenceGenerator.SEgeneratorWithIRoption(RFmatrix1,RFmatrix2,...
            acquiredkspace,acquiredFOV,receiverBW,G_max,TR,TE,dt,gamma,...
            conn_localdb,experiment_id,pulseq_id,partialFourierStruct,0,[],[],[],...
            struct_pulseq);
    end
    info = pulseSequenceGenerator.createInfoStruct_SE(kspace,FOV,dt,...
        receiverBW,TR,TE,kspace_times,f0,sliceThickness,acquiredkspace,...
        acquiredFOV,encodedkSpace);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_SE(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times_new,f0,sliceThickness,...
            acquiredkspace,acquiredFOV,encodedkSpace);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2)),' ',partialFourierText];
    printForImage.thirdLine     = ['TR/TE: ',num2str(TR*1000),'/',num2str(TE*1000),'ms'];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        if strcmp(pulseSeqFamilyName,'IR-SE')
            printForImage.fourthLine = ['TI: ',num2str(str2num(struct_pulseq.ir_time)*1000),'ms',' - SENSE: ',struct_pulseq.rfactor];
        else
            printForImage.fourthLine    = ['SENSE: ',struct_pulseq.rfactor];
        end
    else
        if strcmp(pulseSeqFamilyName,'IR-SE')
            printForImage.fourthLine = ['TI: ',num2str(str2num(struct_pulseq.ir_time)*1000),'ms'];
        else
            printForImage.fourthLine    = '';
        end
    end

elseif strcmp(pulseSeqFamilyName,'TSE') || strcmp(pulseSeqFamilyName,'IR-TSE') || ...
        strcmp(pulseSeqFamilyName,'SS-FSE')

    TR              = str2num(struct_pulseq.tr);
    
    % If SS-FSE, the ETL should be equal to the total number of acquired
    % lines
    if strcmp(pulseSeqFamilyName,'SS-FSE')
        ETL = acquiredkspace(1,2);
    else
        ETL = str2num(struct_pulseq.etl);
    end
    
    sliceThickness  = str2num(struct_pulseq.sliceThickness);
    
    if isfield(struct_pulseq,'fast')
        fast = str2num(struct_pulseq.fast);
    else
        fast = 0;
    end
    
    gamma           = 42.56*10^6;                           % Hz/T
    %dt              = 10^-6;                                % sec
    f0              = 64000000;                             % in Hz

    RFmatrix1        = zeros(6,1);
    RFmatrix1(1,:)   = str2num(...
        struct_pulseq.RFduration);       % RF duration 
    RFmatrix1(2,:)   = 2;                % RF sinc cycles
    RFmatrix1(3,:)   = 90;               % RF rotation angle
    RFmatrix1(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
    RFmatrix1(5,:)   = [1];              % starting timestep of RF
    RFmatrix1(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis

    RFmatrix2        = zeros(6,1);
    RFmatrix2(1,:)   = str2num(...
        struct_pulseq.RFduration);       % RF duration 
    RFmatrix2(2,:)   = 2;                % RF sinc cycles
    if isfield(struct_pulseq,'refoc_pulse_angle')
        RFmatrix2(3,:)   = str2num(...
            struct_pulseq.refoc_pulse_angle);
    else
        RFmatrix2(3,:)   = 180;          % RF rotation angle
    end
    RFmatrix2(4,:)   = sliceThickness;   % desired slice thickness  %0.015(initial)
    RFmatrix2(5,:)   = [1];              % starting timestep of RF
    RFmatrix2(6,:)   = [1];              % 1 is for x-axis, 2 is for y-axis
    
    receiverBW      = str2num(struct_pulseq.BW);            % 200*10^3; % in Hz (this is a feature of the hardware)
    G_max           = 40*10^-3;                             % in T/m

    TR              = str2num(struct_pulseq.tr);
    TE              = str2num(struct_pulseq.te);
    
    slabThickness   = sliceThickness;
    
    [pulse_sequence,isInKspace,pulse_sequence_TR,soft_crushers,N_pulse,N_TR,...
        kspace_times,BW_rf,TR_times,TRmin,orderOfKspacelines,...
        TEeffective,tirl,encodedkSpace] = pulseSequenceGenerator.TSEgenerator(RFmatrix1,RFmatrix2,...
        acquiredkspace,acquiredFOV,receiverBW,G_max,TR,TE,dt,gamma,ETL,conn_localdb,...
        experiment_id,pulseq_id,struct_pulseq,pulseSeqFamilyName,partialFourierStruct);
    info = pulseSequenceGenerator.createInfoStruct_TSE(kspace,FOV,dt,receiverBW,TR,TE,...
        kspace_times,f0,sliceThickness,orderOfKspacelines,acquiredkspace,...
        acquiredFOV,encodedkSpace);
    
    %% ADD FAST ALGORITHM
    if fast == 1
        [pulse_sequence,kspace_times_new,isInKspace,soft_crushers] = ...
            pulseSequenceGenerator.pulseOptimization(pulse_sequence,N_pulse,...
            info.pulseSequence.kspace,soft_crushers);
        N_pulse = size(pulse_sequence,2);
        info = pulseSequenceGenerator.createInfoStruct_TSE(kspace,FOV,dt,...
            receiverBW,TR,TE,kspace_times,f0,sliceThickness,...
            orderOfKspacelines,acquiredkspace,acquiredFOV,encodedkSpace);
    end
    printForImage.firstLine     = ['FOV: ',num2str(FOV(1,1)*1000),'x',num2str(FOV(1,2)*1000),'mm'];
    printForImage.secondLine    = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];
    printForImage.thirdLine     = ['TR/TEeff: ',num2str(TR*1000),'/',num2str(TEeffective*1000),'ms - ETL: ',num2str(ETL)];
    if strcmp(struct_pulseq.parallelImaging,'sense')
        if strcmp(pulseSeqFamilyName,'IR-TSE') || ...
            (strcmp(pulseSeqFamilyName,'TSE') && isfield(struct_pulseq,'fatsat') &&...
            strcmp(struct_pulseq.fatsat,'lipidir'))
            printForImage.fourthLine = ['TI: ',num2str(str2num(struct_pulseq.ir_time)*1000),'ms',' - SENSE: ',struct_pulseq.rfactor];
        else
            printForImage.fourthLine = ['SENSE: ',struct_pulseq.rfactor];
        end
    else
        if strcmp(pulseSeqFamilyName,'IR-TSE') || ...
            (strcmp(pulseSeqFamilyName,'TSE') && isfield(struct_pulseq,'fatsat') &&...
            strcmp(struct_pulseq.fatsat,'lipidir'))
            printForImage.fourthLine = ['TI: ',num2str(str2num(struct_pulseq.ir_time)*1000),'ms'];
        else
            printForImage.fourthLine = '';
        end
    end
end
% In GRE-3D, the tirl is the one initially calculated
if ~strcmp(pulseSeqFamilyName,'GRE-3D') 
    if strcmp(pulseSeqFamilyName,'GRE-conc')
        tirl = tirl * str2double(struct_pulseq.NEX);
    else
        tirl = tirl * str2double(struct_pulseq.NEX) * noOfSlices;
    end
end
%% ALARM - CHECK IF SPURIOUS ECHOES WILL APPEAR
% PHASE DIFFERENCE OF CONTIGUOUS ELEMENTS ALARM
% testSpuriousEchoes - define whether Spurious Echoes Check will be 
% performed on the pulse sequence
sqlquery        = ['SELECT selected_value FROM ',...
    'edt_tool_local.global_configuration WHERE',...
    ' name=''testSpuriousEchoes'''];
sqlquery_results        = exec(conn_localdb, sqlquery);
testSpuriousEchoesInfo  = fetch(sqlquery_results);
testSpuriousEchoes      = str2double(testSpuriousEchoesInfo.Data{1,1});
    
elementsize = [0.001,0.001,0.001];    
[elementsize_new,~,alarm_count] = ...
    pulseSequenceGenerator.testAccumPhase(...
    TRmin,dt,pulse_sequence_TR(4:6,:),...
    pulse_sequence_TR(7:8,:),elementsize,gamma);
if alarm_count>0
    if advancedNotifications
        msg_spEchoes = ['Artificial artifacts (related to the simulation ',...
            'model) may appear and deteriorate the quality of the image. If ',...
            'yes, please try to deactivate the menu item Imaging -> Advanced -> ',...
            'Simulate slice thickness (if active) or decrease the gradient strength.'];
        eduTool.frontend.errorAndDBentry(conn_localdb,msg_spEchoes,'info',...
            experiment_id,pulseq_id)
    end
end
disp(['The alarm was enabled ',num2str(alarm_count),' times.'])
disp(['The elementsize should be ',num2str(elementsize_new(1,1)),' x ',...
    num2str(elementsize_new(1,2)),' x ',num2str(elementsize_new(1,3)),...
    ' meters (m)'])
if (isnan(testSpuriousEchoes) || testSpuriousEchoes) && alarm_count> 0
    errorFlag = 0;
    msg = ['WARNING: The current configuration of the pulse ',...
            'sequence will impose the creation of spurious echoes ',...
            'due to the discretized nature of the anatomical model. ',...
            'Please perform the following actions:<br/><br/>'];            
    if elementsize_new(1,1)<elementsize(1,1)
        msg = [msg,'Consider decreasing the strength of the ',...
            'gradients on the readout direction.<br/>'];
        errorFlag = 1;
    end
    if elementsize_new(1,2)<elementsize(1,2)
        msg = [msg,'Consider decreasing the strength of the gradients ',...
            'on the phase encoding direction.<br/>'];
        errorFlag = 1;
    end
    if elementsize_new(1,3)<elementsize(1,3) &&...
            performance_gridZ_sliceThickness==1
        % performance_gridZ_sliceThickness==1  means that the 
        % anatomical model consists of a thick voulme of spins
        msg = [msg,'Consider decreasing the strength of the gradients ',...
            'on the slice selection direction.<br/>'];
        errorFlag = 1;
    end
    if errorFlag
        status = 'cancelled-error';
        msg = [msg,'OR consider deactivating the Simulation -> Advanced -> Check for spurious echoes test.<br/>'];
        eduTool.frontend.errorAndDBentry(conn_localdb,msg,status,experiment_id,pulseq_id)
    end
end