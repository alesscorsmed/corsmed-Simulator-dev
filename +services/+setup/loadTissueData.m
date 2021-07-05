function [tissueData] = loadTissueData(experimentData)
%
% SERVICES.SETUP.LOADTISSUEDATA
%
%	Function that loads tissue data from (REDIS) json
%
% INPUT
%   experimentData   structure with json loaded parameters from Redis
%
% OUTPUT
%   tissueData   update tissue data with new values
%
%========================  CORSMED AB © 2020 ==============================
%

numProperties   = 7;
numTissues      = max([experimentData.tissues.tissueId]);
tissueValues    = zeros(numTissues,numProperties);
for ii = length(experimentData.tissues):-1:1
    tid = experimentData.tissues(ii).tissueId;
    % get the new Tissue properties: convert from ms to s in T1/T2
    tissueNames{tid}    = experimentData.tissues(ii).name;
    tissueValues(tid,1) = experimentData.tissues(ii).t1*1e-3;
    tissueValues(tid,2) = experimentData.tissues(ii).t2*1e-3;
    tissueValues(tid,3) = experimentData.tissues(ii).pd;
    tissueValues(tid,4) = experimentData.tissues(ii).cs;
    tissueValues(tid,5) = experimentData.tissues(ii).density;
    tissueValues(tid,6) = experimentData.tissues(ii).electricSigma;
    tissueValues(tid,7) = experimentData.tissues(ii).susc;
end
% set PD of zero T1/T2 to zero
tissueValues((tissueValues(:,1) < 1e-4)|(tissueValues(:,2) < 1e-4), 3) = 0;

% assign
tissueData.tissueValues = tissueValues;
tissueData.tissueNames  = tissueNames;
