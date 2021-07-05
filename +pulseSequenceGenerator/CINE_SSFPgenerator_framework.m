function [pulse_sequence_CINE,soft_crushers_CINE,isInKspace_CINE,...
    N_pulse_CINE,info] = CINE_SSFPgenerator_framework(matfile,segments,cardiacCycleDur,...
    acquireAllphases,selectedPhases,plotON)
% This function creates a CINE SSFP pulse sequence based on a bSSFP pulse 
% sequence previously developed. The user can define a steady cardiac cycle
% duration, the number of segments per kspace and it can select specific
% cardiac phases to be acquired. The design of this pulse sequence 
% maintains the steady state through continuous pulsing. Moreover, to 
% achieve and maintain the steady state, RF and gradient pulsing are 
% applied for one preparatory heartbeat.
%
% INPUTS:
% - matfile: the path of the .mat file that keeps the bSSFP pulse sequence
% - segments: the number of segments the kspace is divided to
% - cardiacCycleDur: the duration of the cardiac cycle in sec. This
% duration is retained throughout the entire pulse sequence
% - acquireAllphases: if 1, all cardiac phases will be acquired and the 
% selectedPhases variable will be ignored. If 0, only the cardiac phases
% defined in the selectedPhases variable will be acquired. For these
% phases, the isInKspace will be filled and the kspace matrix will be
% formed.
% - selectedPhases: the cardiac phases that will be acquired
% - plotON: if 1, the figure of the pulse sequence will be plotted
%
% OUTPUTS: 
% - pulse_sequence: the regular pulse_sequence matrix
% - soft_crushers: the regular soft_crushers matrix
% - isInKspace: the regular isInKspace matrix
% - N_pulse: the regular N_pulse matrix
% - info: the regular info matrix. It also contains the info.pulseSequence.
% TR_ALL_times that holds the times of all the TRs, including those for
% which the receiver is not ON. It also contains the info.pulseSequence.
% TR_ALL_acquisition that shows for which TRs the receiver is ON (1) or OFF
% (0)
%
% EXAMPLE:
% [pulse_sequence_CINE,soft_crushers_CINE,isInKspace_CINE,N_pulse_CINE,info] = ...
% pulseSequenceGenerator.CINE_SSFPgenerator_framework('bSSFP_128x128.mat',...
% 16,1.2,0,[3,18],1);

%% Create or load bSSFP
tic
load(matfile,'pulse_sequence','soft_crushers','isInKspace','N_pulse',...
    'info','dt');

%% Remove the a/2 block which has a duration of TE
a2_block    = pulse_sequence(:,1:info.pulseSequence.TR_times(1,1)-1);

timeWithinTRTillACQ_tmstps  = info.pulseSequence.kspace(1,2)-...
    info.pulseSequence.TR_times(1,1);
ACQ_object_duration         = info.pulseSequence.kspace(1,3);

%% Check if the selectedPhases belong to the range of available phases
kspace              = [info.pulseSequence.Nx,info.pulseSequence.Ny];
dt                  = info.pulseSequence.dt;
TE                  = info.pulseSequence.EchoTime;
TR                  = info.pulseSequence.RepetitionTime;
TR_times            = info.pulseSequence.TR_times;
viewsPerSegment     = floor(kspace(1,2)/segments);
cardiacPhases       = floor((cardiacCycleDur-size(a2_block,2)*dt)/...
    (viewsPerSegment*TR));

if ~acquireAllphases
    
    if max(selectedPhases)>cardiacPhases
        
        msg = ['You have requested to acquire a phase that it is not possible. ',...
            'The maximum available phase is ',...
            num2str(cardiacPhases),'. Please consider increasing the cardiacCycleDur, ',...
            ' increasing segments or use fewer selectedPhases.'];
        
        error(msg)
        
    end
    
end

%% Build preparatory heartbeat

totalCardiacCycles  = 1 + ceil(kspace(1,2)/viewsPerSegment); % +1 to 
% achieve and maintain the steady state, RF and gradient pulsing are 
% applied for one preparatory heartbeat

% Keep the TR that will be copied several times
toBeCopiedTR                    = pulse_sequence(:,TR_times(1,1):TR_times(1,2));
toBeCopiedTR_isInKspace         = isInKspace(1,TR_times(1,1):TR_times(1,2));
toBeCopiedTR_isInKspace_inds    = find(toBeCopiedTR_isInKspace);
toBeCopiedTR_isInKspace         = toBeCopiedTR_isInKspace(1,...
    toBeCopiedTR_isInKspace_inds);
tmstpsTillACQwithinTR           = info.pulseSequence.kspace(1,2) - TR_times(1,1);
tmstpsACQwithinTR               = info.pulseSequence.kspace(1,3);

% Find the RF phase of the first RF pulse (not the a/2 pre-pulse)
absoluteMinRFphaseValue = abs(min(toBeCopiedTR(2,:)));
absoluteMaxRFphaseValue = abs(max(toBeCopiedTR(2,:)));
if absoluteMaxRFphaseValue>absoluteMinRFphaseValue
    initialRFphase = absoluteMaxRFphaseValue;
else
    initialRFphase = -absoluteMinRFphaseValue;
end

% Find the time when the TR for the acquisition of the first kspace line
% will start
timeBeforeStartAcq_tmsps    = floor(cardiacCycleDur/dt);

% Find how many times the toBeCopiedTR should be copied
noCopiesTR = ceil((timeBeforeStartAcq_tmsps - size(a2_block,2))/round(TR/dt));

timeBeforeStartAcq          = size(a2_block,2)*dt + noCopiesTR*TR;
timeBeforeStartAcq_tmstps   = size(a2_block,2) + noCopiesTR*round(TR/dt);

% TR_all_times_CINE holds the start and end times of all the TRs (not only
% the ones where acquisition occurs. The third column holds 0s or 1s that
% notifies when this TR is used for acquisition
TR_all_times_CINE   = [];

TR_times_start      = cumsum([1,repmat(round(TR/dt),1,noCopiesTR-1)]);
TR_times_end        = cumsum([round(TR/dt),repmat(round(TR/dt),1,noCopiesTR-1)]);
TR_times_temp       = [TR_times_start',TR_times_end'];

temp_time           = size(a2_block,2);
TR_all_times_CINE   = [TR_all_times_CINE;TR_times_temp+repmat(temp_time,...
    size(TR_times_temp,1),size(TR_times_temp,2))];

TR_all_times_CINE(:,3) = zeros(size(TR_all_times_CINE,1),1);

pulse_sequence_CINE = [a2_block,repmat(toBeCopiedTR,1,noCopiesTR)];
soft_crushers_CINE  = zeros(1,size(pulse_sequence_CINE,2));

TR_times_CINE       = [];
TR_id_CINE          = [];
TR_phase_CINE       = [];
TR_segment_CINE     = [];

TR_temp_index       = 0;

for n = 1:segments

    startTimeToBeCopied = TR_times(((n-1)*viewsPerSegment)+1,1);
    
    if n == segments
        endTimeToBeCopied = TR_times(end,2);
        noTRs = kspace(1,2) - (n-1)*viewsPerSegment;
    else
        endTimeToBeCopied = TR_times(n*viewsPerSegment,2);
        noTRs = viewsPerSegment;
    end
    
    TR_segment_CINE = [TR_segment_CINE;n*ones(noTRs*cardiacPhases,1)];
    
    partToBeCopied_pulse_sequence   = ...
        pulse_sequence(:,startTimeToBeCopied:endTimeToBeCopied);
    partToBeCopied_soft_crushers    = ...
        soft_crushers(:,startTimeToBeCopied:endTimeToBeCopied);
    
    pulse_sequence_CINE             = [pulse_sequence_CINE,...
        repmat(partToBeCopied_pulse_sequence,1,cardiacPhases)];
    soft_crushers_CINE              = [soft_crushers_CINE,...
        repmat(partToBeCopied_soft_crushers,1,cardiacPhases)];
    
    timeBeforeStartCardiacPhase_tmstps = cumsum([timeBeforeStartAcq_tmstps,...
        repmat(size(partToBeCopied_pulse_sequence,2),1,cardiacPhases)]);
    
    TR_id_temp      = 1:noTRs;
    TR_id_temp      = TR_temp_index + TR_id_temp;
    TR_id_CINE      = [TR_id_CINE;repmat(TR_id_temp',cardiacPhases,1)];
    
    TR_temp_index   = TR_temp_index + noTRs;
    
    TR_times_start  = cumsum([1,repmat(round(TR/dt),1,noTRs-1)]);
    TR_times_end    = cumsum([round(TR/dt),repmat(round(TR/dt),1,noTRs-1)]);
    TR_times_temp   = [TR_times_start',TR_times_end'];
    
    for i = 1:cardiacPhases
        
        temp_time           = timeBeforeStartCardiacPhase_tmstps(1,i);
        TR_times_CINE       = [TR_times_CINE;TR_times_temp+repmat(temp_time,...
            size(TR_times_temp,1),size(TR_times_temp,2))];
        
        if ~acquireAllphases 
            
            if ismember(i,selectedPhases)
                
                temp_column = ones(size(repmat(...
            temp_time,size(TR_times_temp,1),size(TR_times_temp,2)),1),1);
        
            else
                
                temp_column = zeros(size(repmat(...
            temp_time,size(TR_times_temp,1),size(TR_times_temp,2)),1),1);
                
            end
            
        else
            
            temp_column = ones(size(repmat(...
                temp_time,size(TR_times_temp,1),size(TR_times_temp,2)),1),1);   
            
        end
        
        TR_all_times_CINE   = [TR_all_times_CINE;[TR_times_temp+repmat(...
            temp_time,size(TR_times_temp,1),size(TR_times_temp,2)),...
            temp_column]];
%         disp(['Cardiac phase:',num2str(i),' - Size TR_all_times_CINE: ',...
%             num2str(size(TR_all_times_CINE,1))])
        
        TR_phase_CINE       = [TR_phase_CINE;repmat(i,noTRs,1)];
        
    end
    
    if n<segments
        
        timeBeforeStartAcq = timeBeforeStartAcq + ...
            size(repmat(partToBeCopied_pulse_sequence,1,cardiacPhases),2)*dt;
        
        if timeBeforeStartAcq < (n+1)*cardiacCycleDur
            
            residualTime    = (n+1)*cardiacCycleDur - timeBeforeStartAcq;
            noCopiesTR      = ceil(residualTime/TR);
            
            startTimeToBeCopied = TR_times(n*viewsPerSegment,1);
            endTimeToBeCopied   = TR_times(n*viewsPerSegment,2);            

            partToBeCopied_pulse_sequence   = ...
                pulse_sequence(:,startTimeToBeCopied:endTimeToBeCopied);
            partToBeCopied_soft_crushers    = ...
                zeros(1,size(partToBeCopied_pulse_sequence,2));

            pulse_sequence_CINE             = [pulse_sequence_CINE,...
                repmat(partToBeCopied_pulse_sequence,1,noCopiesTR)];
            soft_crushers_CINE              = [soft_crushers_CINE,...
                repmat(partToBeCopied_soft_crushers,1,noCopiesTR)];
            
            TR_times_start                  = cumsum([1,...
                repmat(round(TR/dt),1,noCopiesTR-1)]);
            TR_times_end                    = cumsum([round(TR/dt),...
                repmat(round(TR/dt),1,noCopiesTR-1)]);
            TR_times_temp                   = [TR_times_start',...
                TR_times_end'];
            
            temp_time           = TR_all_times_CINE(end,2);
            TR_all_times_CINE   = [TR_all_times_CINE;...
                [TR_times_temp+repmat(temp_time,...
                size(TR_times_temp,1),size(TR_times_temp,2)),...
                zeros(size(repmat(temp_time,...
                size(TR_times_temp,1),size(TR_times_temp,2)),1),1)]];
            
            timeBeforeStartAcq  = timeBeforeStartAcq + ...
                size(repmat(partToBeCopied_pulse_sequence,1,noCopiesTR),2)*dt;
            
        end
        
        timeBeforeStartAcq_tmstps = round(timeBeforeStartAcq/dt);
    
    end
    
end

ACQ_CINE_start = TR_times_CINE(:,1) + timeWithinTRTillACQ_tmstps;

% Build the kspace matrix
kspace_CINE             = zeros(size(ACQ_CINE_start,1),16);
kspace_CINE(:,1)        = [1:size(ACQ_CINE_start,1)]';
kspace_CINE(:,2)        = ACQ_CINE_start;
kspace_CINE(:,3)        = ACQ_object_duration*ones(size(ACQ_CINE_start,1),1);
kspace_CINE(:,4)        = zeros(size(ACQ_CINE_start,1),1);
kspace_CINE(:,5)        = ACQ_object_duration*ones(size(ACQ_CINE_start,1),1);
kspace_CINE(:,6)        = TR_id_CINE;
kspace_CINE(:,7:10)     = ones(size(ACQ_CINE_start,1),4);
kspace_CINE(:,11)       = TR_phase_CINE;
kspace_CINE(:,12:13)    = ones(size(ACQ_CINE_start,1),2);
kspace_CINE(:,14)       = TR_segment_CINE;
kspace_CINE(:,15)       = zeros(size(ACQ_CINE_start,1),1);

%%
isInKspace_CINE     = zeros(1,size(pulse_sequence_CINE,2));
index               = 0;
maxValueIsInKspace  = 0;
for iTR = 1:size(TR_all_times_CINE,1)
    
    pulse_sequence_CINE(2,TR_all_times_CINE(iTR,1):TR_all_times_CINE(iTR,2)) = ...
        cosd((iTR-1)*180)*toBeCopiedTR(2,:);
    TR_RFphase = cosd((iTR-1)*180)*initialRFphase;    
    
    if TR_all_times_CINE(iTR,3) == 1
        
        index = index + 1;
        kspace_CINE(index,16) = TR_RFphase;
        isInKspace_CINE(1,TR_all_times_CINE(iTR,1)+tmstpsTillACQwithinTR:...
            TR_all_times_CINE(iTR,1)+tmstpsTillACQwithinTR+tmstpsACQwithinTR-1) = ...
            maxValueIsInKspace + toBeCopiedTR_isInKspace;
        
        maxValueIsInKspace = max(isInKspace_CINE(1,:));
        
    end
    
end

%% Clear kspace_CINE from unnecessary cardiac phases
if ~acquireAllphases 
    
    for i = 1:cardiacPhases
            
        if ~ismember(i,selectedPhases)

            temp_inds = find(kspace_CINE(:,11)==i);
            kspace_CINE(temp_inds,:) = [];

        end

    end
    
    kspace_CINE(:,1) = [1:size(kspace_CINE(:,1),1)]';
        
end

%% Fix the info structure
info.pulseSequence.TR_times             = TR_all_times_CINE(find(...
    TR_all_times_CINE(:,3)),1:2);
info.pulseSequence.TR_ALL_times         = TR_all_times_CINE(:,1:2);
info.pulseSequence.TR_ALL_acquisition   = TR_all_times_CINE(:,3);
info.pulseSequence.kspace               = kspace_CINE;
info.pulseSequence.basebFFSP            = matfile;

N_pulse_CINE = size(pulse_sequence_CINE,2);

toc
%% PLOT
if plotON == 1
    timeaxis = [1:size(pulse_sequence_CINE,2)]*dt;

    a1 = subplot(7,1,1);
    plot(timeaxis,pulse_sequence_CINE(1,:),'r')
    ylabel('RF magnitude') % in Tesla

    a2 = subplot(7,1,2);
    plot(timeaxis,pulse_sequence_CINE(2,:),'r')
    ylabel('RF phase') % in radians

    a3 = subplot(7,1,3);
    plot(timeaxis,pulse_sequence_CINE(3,:),'r')
    ylabel('RF Freq.') % in Hz

    a4 = subplot(7,1,4);
    plot(timeaxis,pulse_sequence_CINE(4,:))
    ylabel('Gx') % in T/m

    a5 = subplot(7,1,5);
    plot(timeaxis,pulse_sequence_CINE(5,:))
    ylabel('Gy') % in T/m

    a6 = subplot(7,1,6);
    plot(timeaxis,pulse_sequence_CINE(6,:))
    ylabel('Gz') % in T/m

    a7 = subplot(7,1,7);
    plot(timeaxis,isInKspace_CINE,'k')
    ylabel('isInKspace')
    set(gcf,'Name','Pulse Sequence')

    linkaxes([a1 a2 a3 a4 a5 a6 a7],'x');
    set(gca,'XTick',[])
    
%     figure()
%     plot(pulse_sequence_CINE(1,:),'r')
%     hold on
%     plot(TR_times_CINE(:,1)',pulse_sequence_CINE(1,TR_times_CINE(:,1)'),'gs')
%     hold on
%     plot(TR_times_CINE(:,2)',pulse_sequence_CINE(1,TR_times_CINE(:,2)'),'gs')
%     ylabel('RF magnitude') % in Tesla
end