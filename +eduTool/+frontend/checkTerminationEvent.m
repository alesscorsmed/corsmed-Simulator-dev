function restartMCR = checkTerminationEvent(connLocalDB)
%
% EDUTOOL.FRONTEND.EXPECTEXPERIMENT
%
%	repetitive function until an experiment is ready to run, returning
%	first info
%
% INPUT
%
%   connlocalDB             the local connection object
%
% OUTPUT
%   expControl   expControl structure
%
%========================  CORSMED AB © 2020 ==============================
%

sqlquery = ['SELECT selected_value FROM edt_tool_local.global_configuration ',...
            'WHERE name=''edutool_version'''];
        
        
sqlqueryResults = exec(connLocalDB, sqlquery);
fetchResults    = fetch(sqlqueryResults);
eduToolVersion  = fetchResults.Data{1,1};

if strcmp(eduToolVersion,'v1')
    restartMCR = 1;
else
    restartMCR = 0;
end