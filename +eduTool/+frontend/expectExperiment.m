function experiment = expectExperiment(connlocalDB)
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

sqlquery = ['SELECT id,remoteDB_id,status,info,reconstructor,recon_info',...
            ',pulseq_mat_file_name,pulseq_id FROM edt_tool_local.experiments ',...
            'WHERE status=''kernel'''];
        
        
sqlquery_results = exec(connlocalDB, sqlquery);
experiment = fetch(sqlquery_results);