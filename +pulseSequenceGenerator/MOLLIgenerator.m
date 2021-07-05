function [pulse_sequence_MOLLI,times_fitting_single_point,tirl] = ...
    MOLLIgenerator(pulse_sequence_READOUT,info_READOUT,...
    pulse_sequence_IR,cardiac_cycle_duration,TD,TI,dt,MOLLI_scheme,pause_cc,...
    TE,TR,RF_duration,acquiredkspace,conn_localdb,experiment_id,pulseq_id)
% TD defines the time till the onset of the bSSFP readout. It actually
% describes the timepoint of the cardiac cycle when the bSSFP readout will
% be utilized. For example, a TD equal to 0.7s would place the readout at
% the end-diastole of a cardiac cycle of duration 1sec

cardiac_cycle_duration_timesteps = round(cardiac_cycle_duration/dt);

TD_timesteps    = round(TD/dt);

TI_sec          = TI/1000;  % in sec
TI_timesteps    = round(TI_sec/dt);

int_pause_initial           = pause_cc*cardiac_cycle_duration;  % pause among the LL experiments (in sec)
int_pause_timesteps_initial = round(int_pause_initial/dt);

MOLLI_total_duration_timesteps = sum(MOLLI_scheme)*cardiac_cycle_duration_timesteps + ...
    (size(MOLLI_scheme,2)-1)*int_pause_timesteps_initial;

pulse_sequence  = zeros(8,MOLLI_total_duration_timesteps);
soft_crushers   = zeros(1,MOLLI_total_duration_timesteps);
isInKspace      = zeros(1,MOLLI_total_duration_timesteps);

receiver_phase      = [];
kspace_timesteps    = [];

isInKspace_offset           = 0;
total_timesteps_till_now    = 0;
index                       = 0;

for scheme = 1:size(MOLLI_scheme,2)
    
    for LL = 1:MOLLI_scheme(1,scheme)
        
        index = index + 1;
        TI_LL(1,index) = TI(1,index) + (LL-1)*cardiac_cycle_duration*1000;  % in msec
        
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

        if LL~=1
            ccduration = TI(1,index) - TI(1,index-1);
            ccduration = ccduration/1000;  % convert it to sec
            ccduration_timesteps = round(ccduration/dt);
            TD_residual_timesteps = ccduration_timesteps - ...
                round(pulse_sequence_READOUT.N_pulse);
            
            total_timesteps_till_now = total_timesteps_till_now + TD_residual_timesteps;
            
        end
        
        start_timestep = total_timesteps_till_now + 1;
        
        if LL == 1
            
            if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
                
                % If the number of PE steps is an even number, it is hard 
                % to define the central kspace line. This is why the 
                % almost_center_bSSFP_readout is utilized. However, this
                % may cause a discrepancy between the prescribed TI and the 
                % actual TI. To remove the discrepancy, the TI is
                % calculated from the end of the IR pulse till the 
                % almost_center_bSSFP_readout and not till the 
                % center_bSSFP_readout (see #Check-here-tag) 
                times_fitting_single_point(1,index) = total_timesteps_till_now + ...
                    round(TD(1,index)/dt) + round(almost_center_bSSFP_readout);
                
            else
                times_fitting_single_point(1,index) = total_timesteps_till_now + ...
                    round(TD(1,index)/dt) + round(center_bSSFP_readout);
            end
            
            % timesteps from the R-wave triggering (zero timepoint) till start
            % of the readout
            timesteps_before_readout = TD_timesteps(1,index);
            
            temp_pulse_sequence = [zeros(8,TD_timesteps(1,index)),...
                pulse_sequence_READOUT.pulse_sequence];
            
            % #Check-here-tag
            if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
                timesteps_before_IR = round(TD(1,index)/dt) + ...
                    round(almost_center_bSSFP_readout) - ...
                    TI_timesteps(1,index)-...
                    size(pulse_sequence_IR.pulse_sequence,2);
            else
                timesteps_before_IR = round(TD(1,index)/dt) + ...
                    round(center_bSSFP_readout) - ...
                    TI_timesteps(1,index)-...
                    size(pulse_sequence_IR.pulse_sequence,2);
            end
            
            temp_pulse_sequence(:,timesteps_before_IR+1:timesteps_before_IR+...
                size(pulse_sequence_IR.pulse_sequence,2)) = ...
                pulse_sequence_IR.pulse_sequence;
             
            temp_soft_crushers = [zeros(1,TD_timesteps(1,index)),...
                pulse_sequence_READOUT.soft_crushers]; 
            
            % Add crushers before and after the IR pulse
            temp_soft_crushers(1,timesteps_before_IR-1) = 1;
            temp_soft_crushers(1,timesteps_before_IR+size(pulse_sequence_IR.pulse_sequence,2)+1) = 1;
            
            temp_isInKspace = [zeros(1,TD_timesteps(1,index)),...
                pulse_sequence_READOUT.isInKspace];
            
        else
            
            if mod(size(pulse_sequence_READOUT.kspace_times,1),2) == 0
                times_fitting_single_point(1,index) = total_timesteps_till_now + ...
                    almost_center_bSSFP_readout;
            else
                times_fitting_single_point(1,index) = total_timesteps_till_now + ...
                    center_bSSFP_readout;
            end
            
            timesteps_before_readout = 0;
            
            temp_pulse_sequence = pulse_sequence_READOUT.pulse_sequence;
            temp_soft_crushers  = pulse_sequence_READOUT.soft_crushers;        
            temp_isInKspace     = pulse_sequence_READOUT.isInKspace;
            
        end        
        
        temp1 = find(temp_isInKspace);        
        temp_isInKspace(temp1) = temp_isInKspace(temp1) + isInKspace_offset;
        isInKspace_offset = max(temp_isInKspace);
        
        temp_kspace_timesteps       = pulse_sequence_READOUT.kspace_times;
        temp_kspace_timesteps(:,1)  = temp_kspace_timesteps(:,1) + total_timesteps_till_now + timesteps_before_readout;
%         temp_kspace_timesteps(:,2)  = temp_kspace_timesteps(:,2) + total_timesteps_till_now + timesteps_before_readout - ...
%             temp_kspace_timesteps(:,1) + 1;
        temp_kspace_timesteps(:,2)  = temp_kspace_timesteps(:,2) + total_timesteps_till_now + timesteps_before_readout;
        
        pulse_sequence(:,start_timestep:start_timestep+size(temp_pulse_sequence,2)-1)   = temp_pulse_sequence;        
        soft_crushers(:,start_timestep:start_timestep+size(temp_soft_crushers,2)-1)     = temp_soft_crushers;
        isInKspace(:,start_timestep:start_timestep+size(temp_isInKspace,2)-1)           = temp_isInKspace;
                
        kspace_timesteps = [kspace_timesteps;temp_kspace_timesteps];
        receiver_phase  = [receiver_phase;info_READOUT.pulseSequence.kspace(:,16)];
        
        total_timesteps_till_now = total_timesteps_till_now + size(temp_pulse_sequence,2);
    
    end

    if scheme ~= size(MOLLI_scheme,2)
        ccduration1 = TI(1,sum(MOLLI_scheme(1,1:scheme))) - ...
            TI(1,sum(MOLLI_scheme(1,1:scheme))-1);
        ccduration2 = TI(1,sum(MOLLI_scheme(1,1:scheme))+2) - ...
            TI(1,sum(MOLLI_scheme(1,1:scheme))+1);
        ccduration = mean(ccduration1,ccduration2);
    else
        ccduration = TI(1,sum(MOLLI_scheme(1,1:scheme))) - ...
            TI(1,sum(MOLLI_scheme(1,1:scheme))-1);
    end
    
    ccduration = ccduration/1000;  % convert it to sec
    ccduration_timesteps = round(ccduration/dt);
    TD_residual_timesteps = ccduration_timesteps - ...
        round(pulse_sequence_READOUT.N_pulse);

    int_pause = pause_cc*ccduration;  % pause among the LL experiments (in sec)

    if scheme ~= size(MOLLI_scheme,2)
        int_pause_timesteps = round(int_pause/dt) + TD_residual_timesteps - round(TD(1,index+1)/dt);
    else
        int_pause_timesteps = round(int_pause/dt) + TD_residual_timesteps;
    end
    
    total_timesteps_till_now = total_timesteps_till_now + int_pause_timesteps;
    
end

% times_fitting has the time duration between the end of the IR pulse and
% the time when the center of the kspace is aqruired (center of bSSFP 
% readout). It comes from the durations the ramp and half the duration of 
% the bSSFP. Previously I was counting the duration of the IR pulse as well
% times_fitting = TI_LL+((N_pulse_ramp+N_pulse_Readout/2)*dt*1000);  % in msec
% disp('NOTE: the time array being used in the MOLLI fitting DOES NOT take into account the duration of the INV pulse')

times_fitting   = TI;  % in msec
N_pulse         = MOLLI_total_duration_timesteps;

time_axis = 0:dt:(MOLLI_total_duration_timesteps-1)*dt;

pulse_sequence_MOLLI.pulse_sequence         = pulse_sequence;
pulse_sequence_MOLLI.pulse_sequence(7,:)    = 1:size(pulse_sequence_MOLLI.pulse_sequence,2);
pulse_sequence_MOLLI.pulse_sequence(8,:)    = 1:size(pulse_sequence_MOLLI.pulse_sequence,2);

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
tirl = cardiac_cycle_duration*(sum(MOLLI_scheme)+(size(MOLLI_scheme,2)-1)*pause_cc);