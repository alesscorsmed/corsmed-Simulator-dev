function errorAndDBentry(conn_localdb,msg,status,experiment_id,pulseq_id)
% status may have one of the following entries:
% cancelled-error: The backend makes some calculations. The configuration 
% that the user has selected fails. The backend generates an error message,
% saves this message in the db and changes the status of the experiment.
% Under this status, the pulse sequence becomes editable again.
% 
% confirm: The backend makes some calculations regarding the 
% experiment. The configuration that the user has selected does not fail 
% but the user should be notified regarding the outcome of the simulation 
% based on the configuration of the experiment (for example, the current 
% configuration may generate spurious echoes and, in order to avoid them, 
% the backend had to increase the resolution of the anatomical model which 
% will, in turn, increase the total simulation time). The backend generates
% a message, saves this message in the db and changes the status of the 
% experiment to "confirm". The user receives a pop-up with the 
% message and two buttons that ask the user to continue or not. If the user
% selects to continue, the simulation runs as usual. If the user selects to
% cancel, the experiment's status changes to cancelled-user.
%
% info: it adds an info point in the notification history. No action is
% needed by the user

if strcmp(status,'confirm') || strcmp(status,'cancelled-error')
    
    fprintf(1, '\n');
    fprintf(1, '\n  %s: %s', status, msg);
    fprintf(1, '\n');

    exec(conn_localdb,['UPDATE edt_tool_local.pulse_sequence SET info=''',...
                msg,''' WHERE exper_id=',num2str(experiment_id)]);
    exec(conn_localdb,['UPDATE edt_tool_local.experiments SET status=''',...
                status,''' WHERE id=',num2str(experiment_id)]);
    if strcmp(status,'confirm')
        msgType = 'warning';
    elseif strcmp(status,'cancelled-error')
        msgType = 'error';
    end
elseif strcmp(status,'info')
    msgType = 'info';
end
sqlquery = ['SELECT seqnum FROM edt_tool_local.pulse_sequence ',...
    'WHERE id=',num2str(pulseq_id)];
sqlquery_results = exec(conn_localdb, sqlquery);
a = fetch(sqlquery_results);
pulseq_seqnum  = a.Data{1,1};
exec(conn_localdb,['INSERT INTO edt_tool_local.notifications_history (description,type) ',...
    'VALUES  (''',['Pulse seq. ',num2str(pulseq_seqnum),': ',msg],''',''',msgType,''')']);
if strcmp(status,'confirm') || strcmp(status,'cancelled-error')
    error('MESSAGE-BY-PLATFORM');
end