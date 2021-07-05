function [instanceID] = updateBackendStartCounts(conn)
%
% EDUTOOL.FRONTEND.UPDATEBACKENDSTARTCOUNTS
%
%	Informs the BackEnd that is starting
%
% INPUT
%
%   conn  The DB connection object
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%

sqlquery = ['SELECT backendStart,instance_id FROM edt_tool_local.instance_info WHERE 1'];
sqlquery_backendStart   = exec(conn, sqlquery);
startIndexResult        = fetch(sqlquery_backendStart);
startIndex              = startIndexResult.Data{1,1};
instanceID              = startIndexResult.Data{1,2};

curs2 = exec(conn, ...
    ['UPDATE edt_tool_local.instance_info SET backendStart=',num2str(startIndex+1),...
    ' WHERE instance_id=''',instanceID,'''']);