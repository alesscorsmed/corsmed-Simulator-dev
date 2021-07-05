function [] = updateExperimentStatus(connLocalDB,experimentID,status)
%
%FRONTEND.UPDATEEXPERIMENTSTATUS
%
% updates the experiment Status in the DB
%
% INPUT
%   status              status we want to change into
%   experimentID        ID of the corresponding experiment
%   connLocalDB         the DB connection object
%
% OUTPUT
%   none
%
% ========================= CORSMED AB @ 2020 ============================
%


curs = exec(connLocalDB, ...
    ['UPDATE experiments SET status=''',status,''' WHERE id=',...
    num2str(experimentID)]);
                        
end