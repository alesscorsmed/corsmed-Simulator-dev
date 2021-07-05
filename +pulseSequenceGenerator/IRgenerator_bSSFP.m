function [pulse_sequence,N_pulseIR,N_IRmodule] = IRgenerator_bSSFP(dt,gamma,...
    TI,IRduration,IRtype,bSSFPCenterOfAcquisition,conn_localdb,...
    experiment_id,pulseq_id)
% N_pulseIR holds the length of the pulse_sequence, after the application
% of the fast algorithm
% N_TI holds the length of the pulse sequence in real life

IRcycles    = 2;
IRangle     = 180;

% STEP 1 - Design RF
if strcmp(IRtype,'sinc')    
    [RF,BW,rf_timesteps] = pulseSequenceGenerator.addRFsinc(IRduration,IRcycles,...
        IRangle,dt,gamma); %#ok
end

% In this type of experiment (IR-bSSFP) the inversion time counts from the
% end of the IR pulse till the middle of the central kspace line of the
% bSSFP acquisition. N_TI holds the total timesteps from the end of the IR
% pulse till the start of the bSSFP readout.
N_TI = floor(TI/dt) - bSSFPCenterOfAcquisition;

if N_TI<=0
    msg = ['The selected Inversion Time (TI) is smaller than the minimum available TI (',...
    num2str(bSSFPCenterOfAcquisition*dt),'sec) for the current configuration. Consider',...
    ' increasing the TI or decreasing the duration of the bSSFP readout.'];
    
    messagesDB.errorAndDBentry(conn_localdb,msg,'cancelled-error',experiment_id,pulseq_id);
end

% N_IRmodule holds the total timesteps of the IR module (IR pulse +
% remaining till the start of the bSSFP readout)
N_IRmodule = rf_timesteps + N_TI;

% N_pulseIR holds the length of the compressed version (fast algorithm)...
% of the IR module
N_pulseIR           = rf_timesteps + 1;

pulse_sequence                      = zeros(8,N_pulseIR);
pulse_sequence(1,1:rf_timesteps)    = RF;
pulse_sequence(7,:)                 = 1:size(pulse_sequence,2);
pulse_sequence(8,:)                 = 1:size(pulse_sequence,2);

% Add one extra point that holds the accumulated effect of the pulse 
% sequence till the end of the IR module
pulse_sequence(8,end) = N_IRmodule;