function updateExecTime(expControl,tTotal)
%update the ExecTime_sec in Experiments Table
%
%========================  CORSMED AB © 2020 ==============================
%

if isfield(expControl, 'connLocalDB') && ~isempty(expControl.connLocalDB)
    exec(expControl.connLocalDB,...
        ['UPDATE edt_tool_local.experiments SET execTime_sec=',...
        num2str(tTotal), ...
        ' WHERE id=', num2str(expControl.experimentID)]);
    
end