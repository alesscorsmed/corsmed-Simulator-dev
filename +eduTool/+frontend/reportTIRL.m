function reportTIRL(conn_localdb,tirl,experiment_id)
% tirl stands for "time in real life"

tirl_minutes = floor(tirl/60);
tirl_seconds = round(rem(tirl,60)*100)/100;
if tirl_minutes==0
    if tirl_seconds<1
        tirl_str = [num2str(tirl_seconds),'s'];
    else
        tirl_str = [num2str(round(tirl_seconds*10)/10),'s'];
    end
else
    tirl_str = [num2str(tirl_minutes),'m:',num2str(tirl_seconds),'s'];
end

exec(conn_localdb,['UPDATE edt_tool_local.pulse_sequence SET irl_info=''',...
            tirl_str,''' WHERE exper_id=',num2str(experiment_id)]);