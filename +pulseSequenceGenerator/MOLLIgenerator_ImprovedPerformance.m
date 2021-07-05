function [pulse_sequence_MOLLI,times_fitting_single_point,tirl] = ...
    MOLLIgenerator_ImprovedPerformance(pulse_sequence_READOUT,info_READOUT,...
    pulse_sequence_IR,cardiacCycleDuration,TD,TI,dt,MOLLI_scheme,pause_cc,...
    TE,TR,RF_duration,acquiredkspace,conn_localdb,experiment_id,pulseq_id)
% TD defines the time till the onset of the bSSFP readout. It actually
% describes the timepoint of the cardiac cycle when the bSSFP readout will
% be utilized. For example, a TD equal to 0.7s would place the readout at
% the end-diastole of a cardiac cycle of duration 1sec

cardiacCycleDurationTimesteps = round(cardiacCycleDuration/dt);

TI_sec          = TI/1000;  % in sec
TI_timesteps    = round(TI_sec/dt);

pulse_sequence  = [];
soft_crushers   = [];
isInKspace      = [];

receiver_phase      = [];
kspace_timesteps    = [];

isInKspace_offset           = 0;
total_timesteps_till_now    = 0;
total_indeces_till_now      = 0;
index                       = 0;

bSSFPkspaceTimesteps        = pulse_sequence_READOUT.kspace_times;

for scheme = 1:size(MOLLI_scheme,2)
    
    for LL = 1:MOLLI_scheme(1,scheme)
        
        index = index + 1;
        TI_LL(1,index) = TI(1,index) + (LL-1)*cardiacCycleDuration*1000;  % in msec
        
        if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
            
            center_kspaceline = round(size(pulse_sequence_READOUT.kspace_times,1)/2)-1;
            center_bSSFP_readout = pulse_sequence_READOUT.TR_times(1,1) + ...
                size(pulse_sequence_READOUT.kspace_times,1)*round(TR/dt)/2;  
            almost_center_bSSFP_readout = pulse_sequence_READOUT.TR_times(1,1) + ...
                center_kspaceline * round(TR/dt) + ...
                round(RF_duration/dt/2) + round(TE/dt);
            
        else
            
            center_kspaceline = round(size(pulse_sequence_READOUT.kspace_times,1)/2)-1;
            center_bSSFP_readout = pulse_sequence_READOUT.TR_times(1,1) + ...
                center_kspaceline*round(TR/dt)+...
                round(RF_duration/dt/2)+round(TE/dt);   
            
        end
        
        if center_bSSFP_readout > TI_timesteps(1,1)
            msg = ['The minimum TI should be at least ',num2str(center_bSSFP_readout*dt),...
                'sec. Please modify the TIs of the selected MOLLI scheme.'];

            messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
        end
        
        if LL == 1
            
%             % @@@ modify times_fitting_single_point
%             if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
%                 
%                 % If the number of PE steps is an even number, it is hard 
%                 % to define the central kspace line. This is why the 
%                 % almost_center_bSSFP_readout is utilized. However, this
%                 % may cause a discrepancy between the prescribed TI and the 
%                 % actual TI. To remove the discrepancy, the TI is
%                 % calculated from the end of the IR pulse till the 
%                 % almost_center_bSSFP_readout and not till the 
%                 % center_bSSFP_readout (see #Check-here-tag) 
%                 times_fitting_single_point(1,index) = total_timesteps_till_now + ...
%                     round(TD(1,index)/dt) + round(almost_center_bSSFP_readout);
%                 
%             else
%                 times_fitting_single_point(1,index) = total_timesteps_till_now + ...
%                     round(TD(1,index)/dt) + round(center_bSSFP_readout);
%             end
            
            % #Check-here-tag
            if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
                timesteps_before_IR = round(TD(1,index)/dt) + ...
                    round(almost_center_bSSFP_readout) - ...
                    TI_timesteps(1,index)-...
                    size(pulse_sequence_IR.pulse_sequence,2);
                timestepsTillCenterbSSFP = round(almost_center_bSSFP_readout);
            else
                timesteps_before_IR = round(TD(1,index)/dt) + ...
                    round(center_bSSFP_readout) - ...
                    TI_timesteps(1,index)-...
                    size(pulse_sequence_IR.pulse_sequence,2);                
                timestepsTillCenterbSSFP = round(center_bSSFP_readout);
            end
            
            % The first 2-columns of zeros cover the time from onset till
            % the start of the IR. The second column is for the software 
            % crusher. The second 2-columns of zeros cover the time
            % between the IR and the bSSFP readout. The first column is for 
            % the software crusher at the end of IR. The third column of
            % zeros cover the time from the end of readout till the end of
            % the cc.
            
            timestepsBetweenIRandbSSFP  = round(TD(1,index)/dt) - ...
                timesteps_before_IR - size(pulse_sequence_IR.pulse_sequence,2);
            
            tempPulseSequence   = zeros(8,...
                size(pulse_sequence_IR.pulse_sequence,2)+...
                size(pulse_sequence_READOUT.pulse_sequence,2)+5);
            tempSoftCrushers    = zeros(1,size(tempPulseSequence,2));
            tempIsInKspace      = zeros(1,size(tempPulseSequence,2));            
            
            tempPulseSequence(1:6,:) = [zeros(6,2),...
                pulse_sequence_IR.pulse_sequence(1:6,:),...
                zeros(6,2),...
                pulse_sequence_READOUT.pulse_sequence(1:6,:),...
                zeros(6,1)];
            
            % Update rows 7 and 8 till IR
            tempPulseSequence(7:8,1:2) = [[total_timesteps_till_now+1,total_timesteps_till_now+timesteps_before_IR];...
                [total_timesteps_till_now+timesteps_before_IR-1,total_timesteps_till_now+timesteps_before_IR]];
            currentTimePoint                    = total_timesteps_till_now+timesteps_before_IR;            
            currentIndex                        = 2;
            tempSoftCrushers(1,2)               = 1;
            
            % Update rows 7 and 8 for IR
            tempPulseSequence(7:8,currentIndex+1:currentIndex+size(pulse_sequence_IR.pulse_sequence,2)) = ...
                repmat(currentTimePoint+1:currentTimePoint+size(pulse_sequence_IR.pulse_sequence,2),2,1);            
            currentTimePoint                    = currentTimePoint+size(pulse_sequence_IR.pulse_sequence,2);
            currentIndex                        = currentIndex+size(pulse_sequence_IR.pulse_sequence,2);
            
            % Update rows 7 and 8 till start bSSFP
            tempPulseSequence(7:8,currentIndex+1:currentIndex+2) = ...
                [[currentTimePoint+1,currentTimePoint+2];...
                [currentTimePoint+1,currentTimePoint+timestepsBetweenIRandbSSFP]];          
            currentTimePoint                    = currentTimePoint+timestepsBetweenIRandbSSFP;
            currentIndex                        = currentIndex+2;
            tempSoftCrushers(1,currentIndex-1)  = 1;
            
            times_fitting_single_point(1,index) = total_indeces_till_now + ...
                currentIndex + timestepsTillCenterbSSFP;
            
            %Update the kspace matrix
            temp_kspace_timesteps(:,1)  = bSSFPkspaceTimesteps(:,1) + size(pulse_sequence,2) + currentIndex;
            temp_kspace_timesteps(:,2)  = bSSFPkspaceTimesteps(:,2) + size(pulse_sequence,2) + currentIndex;
            
            % Update rows 7 and 8 for bSSFP
            tempPulseSequence(7:8,currentIndex+1:currentIndex+size(pulse_sequence_READOUT.pulse_sequence,2)) = ...
                repmat(currentTimePoint+1:currentTimePoint+size(pulse_sequence_READOUT.pulse_sequence,2),2,1);
            tempInds    = find(pulse_sequence_READOUT.isInKspace);
            tempArray   = pulse_sequence_READOUT.isInKspace;
            tempArray(1,tempInds) = tempArray(1,tempInds) + isInKspace_offset;
            tempIsInKspace(1,currentIndex+1:currentIndex+size(pulse_sequence_READOUT.pulse_sequence,2)) = ...
                tempArray;
            isInKspace_offset   = max(tempArray);
            currentTimePoint    = currentTimePoint+size(pulse_sequence_READOUT.pulse_sequence,2);
            currentIndex        = currentIndex+size(pulse_sequence_READOUT.pulse_sequence,2);
            
            timestepsTillEndCc                  = total_timesteps_till_now + cardiacCycleDurationTimesteps - currentTimePoint;
            
            % Update rows 7 and 8 till end of cc
            tempPulseSequence(7:8,currentIndex+1) = ...
                [currentTimePoint+1;currentTimePoint+timestepsTillEndCc];            
            currentTimePoint                    = currentTimePoint+timestepsTillEndCc;
            currentIndex                        = currentIndex+1;
            
            total_timesteps_till_now            = currentTimePoint;
            total_indeces_till_now              = total_indeces_till_now + currentIndex;
            
        else
            
%             % @@@ modify times_fitting_single_point
%             if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
%                 times_fitting_single_point(1,index) = total_timesteps_till_now + ...
%                     almost_center_bSSFP_readout;
%             else
%                 times_fitting_single_point(1,index) = total_timesteps_till_now + ...
%                     center_bSSFP_readout;
%             end
            
            tempPulseSequence   = zeros(8,...
                size(pulse_sequence_READOUT.pulse_sequence,2)+2);            
            tempSoftCrushers    = zeros(1,size(tempPulseSequence,2));
            tempIsInKspace      = zeros(1,size(tempPulseSequence,2));
            
            % The first column of zeros cover the time from onset till 
            % the start of the bSSFP readout. The second column of zeros 
            % cover the time from the end of readout till the end of cc
            tempPulseSequence(1:6,:) = [zeros(6,1),...
                pulse_sequence_READOUT.pulse_sequence(1:6,:),...
                zeros(6,1)];
            
            % Update rows 7 and 8 till bSSFP
            tempPulseSequence(7:8,1)  = [total_timesteps_till_now+1;...
                total_timesteps_till_now+round(TD(1,index)/dt)];
            currentTimePoint            = total_timesteps_till_now+round(TD(1,index)/dt);            
            currentIndex                = 1;
            
            times_fitting_single_point(1,index) = total_indeces_till_now + ...
                currentIndex + timestepsTillCenterbSSFP;
            
            % Update the kspace matrix
            temp_kspace_timesteps(:,1)  = bSSFPkspaceTimesteps(:,1) + size(pulse_sequence,2) + currentIndex;
            temp_kspace_timesteps(:,2)  = bSSFPkspaceTimesteps(:,2) + size(pulse_sequence,2) + currentIndex;
            
            % Update rows 7 and 8 for bSSFP
            tempPulseSequence(7:8,currentIndex+1:currentIndex+size(pulse_sequence_READOUT.pulse_sequence,2)) = ...
                repmat(currentTimePoint+1:currentTimePoint+size(pulse_sequence_READOUT.pulse_sequence,2),2,1);
            tempInds    = find(pulse_sequence_READOUT.isInKspace);
            tempArray   = pulse_sequence_READOUT.isInKspace;
            tempArray(1,tempInds) = tempArray(1,tempInds) + isInKspace_offset;
            tempIsInKspace(1,currentIndex+1:currentIndex+size(pulse_sequence_READOUT.pulse_sequence,2)) = ...
                tempArray;
            isInKspace_offset   = max(tempArray);
            currentTimePoint    = currentTimePoint+size(pulse_sequence_READOUT.pulse_sequence,2);
            currentIndex        = currentIndex+size(pulse_sequence_READOUT.pulse_sequence,2);
            
            timestepsTillEndCc                  = total_timesteps_till_now + cardiacCycleDurationTimesteps - currentTimePoint;
            
            % Update rows 7 and 8 till end cc
            tempPulseSequence(7:8,currentIndex+1) = ...
                [currentTimePoint+1;currentTimePoint+timestepsTillEndCc];            
            currentTimePoint                    = currentTimePoint+timestepsTillEndCc;
            currentIndex                        = currentIndex+1;
            
            total_timesteps_till_now            = currentTimePoint;
            total_indeces_till_now              = total_indeces_till_now + currentIndex;
            
        end
        
        pulse_sequence      = [pulse_sequence,tempPulseSequence];
        soft_crushers       = [soft_crushers,tempSoftCrushers];
        isInKspace          = [isInKspace,tempIsInKspace];        
                
        kspace_timesteps    = [kspace_timesteps;temp_kspace_timesteps];
        receiver_phase      = [receiver_phase;info_READOUT.pulseSequence.kspace(:,16)];
    
    end

    if scheme ~= size(MOLLI_scheme,2)
        int_pause = pause_cc(1,scheme)*cardiacCycleDurationTimesteps;  % pause among the LL experiments (in sec)
        
        pulse_sequence              = [pulse_sequence,...
            [zeros(6,1);total_timesteps_till_now+1;total_timesteps_till_now+int_pause]];
        soft_crushers               = [soft_crushers,0];
        isInKspace                  = [isInKspace,0];
        
        total_timesteps_till_now    = total_timesteps_till_now+int_pause;
        total_indeces_till_now      = total_indeces_till_now + 1;
    end
    
end

% times_fitting has the time duration between the end of the IR pulse and
% the time when the center of the kspace is aqruired (center of bSSFP 
% readout). It comes from the durations the ramp and half the duration of 
% the bSSFP. Previously I was counting the duration of the IR pulse as well
% times_fitting = TI_LL+((N_pulse_ramp+N_pulse_Readout/2)*dt*1000);  % in msec
% disp('NOTE: the time array being used in the MOLLI fitting DOES NOT take into account the duration of the INV pulse')

times_fitting   = TI;  % in msec
N_pulse         = size(pulse_sequence,2);

pulse_sequence_MOLLI.pulse_sequence     = pulse_sequence;

pulse_sequence_MOLLI.isInKspace         = isInKspace;
pulse_sequence_MOLLI.soft_crushers      = soft_crushers;
pulse_sequence_MOLLI.N_pulse            = N_pulse;
pulse_sequence_MOLLI.kspace_times       = kspace_timesteps;
pulse_sequence_MOLLI.receiver_phase     = receiver_phase;

pulse_sequence_MOLLI.info               = info_READOUT;

pulse_sequence_MOLLI.info.pulseSequence.kspace = ...
    repmat(pulse_sequence_MOLLI.info.pulseSequence.kspace,sum(MOLLI_scheme),1);
pulse_sequence_MOLLI.info.pulseSequence.kspace(:,2:3)   = [kspace_timesteps(:,1),kspace_timesteps(:,2)-kspace_timesteps(:,1)+1];
pulse_sequence_MOLLI.info.pulseSequence.kspace(:,16)    = receiver_phase;
pulse_sequence_MOLLI.info.pulseSequence.kspace(:,1)     = [1:size(pulse_sequence_MOLLI.info.pulseSequence.kspace,1)]';

temp = repmat(1:sum(MOLLI_scheme),acquiredkspace(1,2),1);
pulse_sequence_MOLLI.info.pulseSequence.kspace(:,10)    = temp(:);

%% Calculate TIRL
tirl = cardiacCycleDuration*(sum(MOLLI_scheme) + sum(pause_cc));