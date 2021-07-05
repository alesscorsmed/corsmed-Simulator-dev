function progressUpdate(expControl)
%
% update progress bar in EduTool
%
%========================  CORSMED AB © 2020 ==============================
%
if isfield(expControl, 'connLocalDB') && ~isempty(expControl.connLocalDB)
    exec(expControl.connLocalDB,...
        ['UPDATE edt_tool_local.experiments_progress SET progress=',...
        num2str(round(expControl.progress)), ...
        ' WHERE exper_id=', num2str(expControl.experimentID)]);
else
    
end
  
%fprintf(1, '\nPROGRESS %d%%\n', round(expControl.progress));