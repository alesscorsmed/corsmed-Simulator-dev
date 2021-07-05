function UIindicators = reportSARtIRL(sarEST,timeIRL,expControl)
%
% update EduTool with SAR and TIRL
%
%========================  CORSMED AB © 2020 ==============================
%

% convert to strings
if (sarEST <=0) || isnan(sarEST)
    strShortSAR = 'N/A';
    strLongSAR  = 'N/A';
else
    % compute average 1 min SAR (sarEst is 1s average)
    strShortSAR = sprintf('%.2f',sarEST*60); % W/Kg/min
    strLongSAR  = sprintf('%.3f',sarEST*60); % W/Kg/min
end
% convert time
timeIRL = seconds(ceil(timeIRL));
timeIRL.Format = 'mm:ss';
strTimeIRL = sprintf('%ss', strrep(sprintf('%s',timeIRL),':', 'm:'));

UIindicators.strTimeIRL     = strTimeIRL;
UIindicators.strShortSAR    = strShortSAR;
UIindicators.strLongSAR     = strLongSAR;

if isfield(expControl, 'connLocalDB') && ~isempty(expControl.connLocalDB)
    % report TIRL
    exec(expControl.connLocalDB,...
        ['UPDATE edt_tool_local.pulse_sequence SET irl_info=''',...
        strTimeIRL,''' WHERE exper_id=',num2str(expControl.experimentID)]);
    % repott SAR
    exec(expControl.connLocalDB,...
        ['UPDATE edt_tool_local.experiments SET sar=''',...
        strShortSAR,''', sar_backend=''',strLongSAR,...
        ''' WHERE id=',num2str(expControl.experimentID)]);
    
end